#include <array>
#include <cmath>
#include <iostream>
#include <set>
#include <type_traits>
#include <vector>

#include "Core/gs_parallel_api.h"
#include "ModelGrid/ModelGrid.h"
#include "ModelGrid/ModelGridEditor.h"
#include "ModelGrid/ModelGridMesher.h"
#include "ModelGrid/ModelGridMeshCache.h"
#include "ModelGrid/ModelGridCell.h"
#include "Math/GSTransformList.h"

#include "src/DenseMeshAPI.h"
#include "src/DummyParallelAPI.h"
#include "src/MeshCode.h"
#include "src/DenseMeshAABBTree.h"
#include "src/GeometryUtils.h"
#include "MeshIO/OBJReader.h"
#include "MeshIO/OBJFormatData.h"
#include "src/BoxGenerator.h"
#include "src/SphereGenerator.h"

// 16x16x16 inside/outside bitmap (4096 bits = 64 int64s)
// Grid dimension matches ModelGridCellData_StandardRST::MaxDimension+1
static constexpr int BitGridN = 16;
static constexpr int BitGridNumBits = BitGridN * BitGridN * BitGridN;
static constexpr int BitGridNumWords = BitGridNumBits / 64;

struct BitGrid
{
    std::array<int64_t, BitGridNumWords> Data = {};

    void SetBit(int Index) {
        Data[Index / 64] |= (int64_t(1) << (Index % 64));
    }

    bool GetBit(int Index) const {
        return (Data[Index / 64] >> (Index % 64)) & 1;
    }

    // Count bits that are the same in both grids
    int CountMatchingBits(const BitGrid& Other) const {
        int Count = 0;
        for (int i = 0; i < BitGridNumWords; ++i) {
            int64_t XorBits = Data[i] ^ Other.Data[i];
            Count += __builtin_popcountll((uint64_t)XorBits);
        }
        return BitGridNumBits - Count;
    }

    // Count bits that differ between the two grids
    int CountDifferentBits(const BitGrid& Other) const {
        int Count = 0;
        for (int i = 0; i < BitGridNumWords; ++i) {
            int64_t XorBits = Data[i] ^ Other.Data[i];
            Count += __builtin_popcountll((uint64_t)XorBits);
        }
        return Count;
    }

    int CountSetBits() const {
        int Count = 0;
        for (int i = 0; i < BitGridNumWords; ++i)
            Count += __builtin_popcountll((uint64_t)Data[i]);
        return Count;
    }

    int CountUnsetBits() const {
        return BitGridNumBits - CountSetBits();
    }

    bool operator==(const BitGrid& Other) const { return Data == Other.Data; }
    bool operator<(const BitGrid& Other) const { return Data < Other.Data; }
};


// Sample a 16x16x16 grid of points inside Bounds at dimension-step midpoints,
// evaluate InsideFunc at each point, and return a BitGrid with results.
BitGrid SampleInsideBitGrid(
    const GS::AxisBox3d& Bounds,
    GS::FunctionRef<bool(const GS::Vector3d&)> InsideFunc)
{
    BitGrid Result;

    GS::Vector3d Min = Bounds.Min;
    GS::Vector3d Size = Bounds.Max - Bounds.Min;

    // Sample at midpoints of the 16 dimension intervals: (i + 0.5) / 16
    int BitIndex = 0;
    for (int iz = 0; iz < BitGridN; ++iz) {
        double tz = ((double)iz + 0.5) / (double)BitGridN;
        for (int iy = 0; iy < BitGridN; ++iy) {
            double ty = ((double)iy + 0.5) / (double)BitGridN;
            for (int ix = 0; ix < BitGridN; ++ix) {
                double tx = ((double)ix + 0.5) / (double)BitGridN;
                GS::Vector3d Point(
                    Min.X + tx * Size.X,
                    Min.Y + ty * Size.Y,
                    Min.Z + tz * Size.Z);
                if (InsideFunc(Point))
                    Result.SetBit(BitIndex);
                BitIndex++;
            }
        }
    }

    return Result;
}


template<typename CellType>
void AppendCellMesh(GS::ModelGridMesher& Mesher, const GS::AxisBox3d& CellBounds,
    const GS::ModelGridMesher::CellMaterials& Materials, GS::DenseMeshBuilder& Builder,
    GS::TransformListd& Transforms, GS::ModelGridMesher::AppendCache& Cache)
{
    if constexpr (std::is_same_v<CellType, GS::MGCell_Slab>)
        Mesher.AppendBox(CellBounds, Materials, Builder, Transforms, Cache);
    else if constexpr (std::is_same_v<CellType, GS::MGCell_Ramp>)
        Mesher.AppendRamp(CellBounds, Materials, Builder, Transforms, Cache);
    else if constexpr (std::is_same_v<CellType, GS::MGCell_Corner>)
        Mesher.AppendCorner(CellBounds, Materials, Builder, Transforms, Cache);
    else if constexpr (std::is_same_v<CellType, GS::MGCell_Pyramid>)
        Mesher.AppendPyramid(CellBounds, Materials, Builder, Transforms, Cache);
    else if constexpr (std::is_same_v<CellType, GS::MGCell_Peak>)
        Mesher.AppendPeak(CellBounds, Materials, Builder, Transforms, Cache);
    else if constexpr (std::is_same_v<CellType, GS::MGCell_Cylinder>)
        Mesher.AppendCylinder(CellBounds, Materials, Builder, Transforms, Cache);
    else if constexpr (std::is_same_v<CellType, GS::MGCell_CutCorner>)
        Mesher.AppendCutCorner(CellBounds, Materials, Builder, Transforms, Cache);
}


// Generate a DenseMesh for a ModelGrid cell with the given typed cell params.
template<typename CellType>
GS::DenseMesh GenerateCellMesh(
    const GS::Vector3d& CellDims,
    const CellType& Cell)
{
    GS::ModelGridMesher Mesher;
    Mesher.Initialize(CellDims);
    GS::ModelGridMesher::AppendCache Cache;
    Mesher.InitAppendCache(Cache);
    GS::DenseMeshBuilder Builder;

    GS::TransformListd Transforms;
    GS::GetUnitCellTransform(Cell, CellDims, Transforms);

    GS::AxisBox3d CellBounds(GS::Vector3d::Zero(), CellDims);
    GS::ModelGridMesher::CellMaterials Materials;
    AppendCellMesh<CellType>(Mesher, CellBounds, Materials, Builder, Transforms, Cache);

    return Builder.ToDenseMesh();
}


struct CellMatchResult
{
    GS::ModelGridCell Cell;
    int MatchingBits = 0;
    int DifferentBits = 0;
};



struct PrecomputedCell
{
    BitGrid BG;
    GS::ModelGridCell Cell;
};


std::vector<PrecomputedCell> PrecomputeCellTable(const GS::Vector3d& CellDims)
{
    GS::AxisBox3d LocalBounds(GS::Vector3d::Zero(), CellDims);
    std::vector<PrecomputedCell> Table;

    // Empty cell
    {
        PrecomputedCell Entry;
        Entry.BG = {};
        Entry.Cell = GS::ModelGridCell::EmptyCell();
        Table.push_back(Entry);
    }

    // Solid cell
    {
        PrecomputedCell Entry;
        for (auto& W : Entry.BG.Data) W = ~int64_t(0);
        Entry.Cell = GS::ModelGridCell::SolidCell();
        Table.push_back(Entry);
    }

    //const unsigned int TryDimensionVals[] = { 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
    const unsigned int TryDimensionVals[] = { 3, 5, 7, 9, 11, 13, 15 };
    //const unsigned int TryDimensionVals[] = { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15 };


    std::set<BitGrid> SeenBitGrids;
    int DuplicateCount = 0;

    auto AddType = [&](auto CellTypeTag) {
        using CellType = decltype(CellTypeTag);
        for (unsigned int Dir = 0; Dir < 6; ++Dir)
            for (unsigned int Rot = 0; Rot < 4; ++Rot)
                for (unsigned int DimX : TryDimensionVals)
                    for (unsigned int DimY : TryDimensionVals)
                        for (unsigned int DimZ : TryDimensionVals) {
                            CellType TypedCell = CellType::GetDefaultCellParams();
                            TypedCell.Params.AxisDirection = Dir;
                            TypedCell.Params.AxisRotation = Rot;
                            TypedCell.Params.DimensionX = DimX;
                            TypedCell.Params.DimensionY = DimY;
                            TypedCell.Params.DimensionZ = DimZ;
                            GS::DenseMesh Mesh = GenerateCellMesh<CellType>(CellDims, TypedCell);
                            auto BG = SampleInsideBitGrid(LocalBounds,
                                [&](const GS::Vector3d& P) { return GS::WindingNumber(Mesh, P) < -0.5; });
                            if (!SeenBitGrids.insert(BG).second) {
                                DuplicateCount++;
                                continue;
                            }
                            PrecomputedCell Entry;
                            Entry.BG = BG;
                            GS::UpdateGridCellFromSubCell(Entry.Cell, TypedCell);
                            Table.push_back(Entry);
                        }
    };

    AddType(GS::MGCell_Slab{});
    AddType(GS::MGCell_Ramp{});
    //AddType(GS::MGCell_Corner{});
    //AddType(GS::MGCell_Pyramid{});
    //AddType(GS::MGCell_Peak{});
    //AddType(GS::MGCell_Cylinder{});
    //AddType(GS::MGCell_CutCorner{});

    std::cout << "  Filtered " << DuplicateCount << " duplicate BitGrids" << std::endl;
    return Table;
}


CellMatchResult FindBestCellMatch_Precomputed(
    const BitGrid& Target,
    const std::vector<PrecomputedCell>& Table)
{
    CellMatchResult Best;

    for (const auto& Entry : Table) {
        int Matching = Entry.BG.CountMatchingBits(Target);
        if (Matching > Best.MatchingBits) {
            Best.Cell = Entry.Cell;
            Best.MatchingBits = Matching;
            Best.DifferentBits = BitGridNumBits - Matching;
        }
        if (Best.MatchingBits == BitGridNumBits)
            return Best;
    }

    return Best;
}


int main()
{
    GS::Parallel::RegisterAPI(GSMakeUniquePtr<GS::DummyParallelAPI>());

    GS::DenseMesh SourceMesh;

    //std::string InputFilePath = "../input/stanford-bunny.obj";
    std::string InputFilePath = "../input/Bauhaus_main.obj";

    // Load source mesh
    // GS::OBJFormatData OBJData;
    // if (!GS::OBJReader::ReadOBJ(InputFilePath, OBJData)) {
    //     std::cerr << "Failed to read " << InputFilePath << std::endl;
    //     return 1;
    // }
    // GS::OBJFormatDataToDenseMesh(OBJData, SourceMesh);

    GS::DenseMeshBuilder Builder;
    GS::BoxGenerator BoxGen;
    BoxGen.Center = GS::Vector3d(2.0, 5.0, 7.0);
    BoxGen.Dimensions = GS::Vector3d(100.0, 100.0, 100.0);
    BoxGen.Generate(Builder);
    SourceMesh = Builder.ToDenseMesh();

    // GS::DenseMeshBuilder Builder;
    // GS::SphereGenerator SphereGen;
    // SphereGen.Radius = 50.0;
    // SphereGen.Slices = 32;
    // SphereGen.Stacks = 24;
    // SphereGen.Generate(Builder);
    // SourceMesh = Builder.ToDenseMesh();


    std::cout << "Source mesh: " << SourceMesh.GetVertexCount() << " vertices, "
              << SourceMesh.GetTriangleCount() << " triangles" << std::endl;

    // Compute cell size from mesh bounds
    constexpr int NumCellsPerDimension = 10;
    GS::AxisBox3d MeshBounds = GS::ComputeBounds(SourceMesh);
    double MaxDim = std::max({MeshBounds.DimensionX(), MeshBounds.DimensionY(), MeshBounds.DimensionZ()});
    double CellSize = MaxDim / NumCellsPerDimension;
    GS::ModelGrid Grid;
    Grid.Initialize(GS::Vector3d(CellSize, CellSize, CellSize));

    GS::ModelGridEditor Editor(Grid);

    GS::Vector3d CellDims(CellSize, CellSize, CellSize);
    std::cout << "Precomputing cell table..." << std::endl;
    auto CellTable = PrecomputeCellTable(CellDims);
    std::cout << "Precomputed " << CellTable.size() << " cell entries ("
              << (CellTable.size() * sizeof(CellTable[0]) + 512) / 1024 << " KB)" << std::endl;

    GS::DenseMeshAABBTree SourceTree;
    SourceTree.Build(SourceMesh);

    auto MeshInside = [&](const GS::Vector3d& P) {
        return SourceTree.FastWindingNumber(P) > 0.5;
    };

    // Compute grid index range covering the mesh bounds
    bool bIsInGrid = false;
    GS::Vector3i MinCell = Grid.GetCellAtPosition(MeshBounds.Min, bIsInGrid);
    GS::Vector3i MaxCell = Grid.GetCellAtPosition(MeshBounds.Max, bIsInGrid);

    constexpr int SkipBitThreshold = 3;

    GS::Vector3i CellCounts = MaxCell - MinCell + GS::Vector3i(1, 1, 1);
    std::cout << "Cell range: " << CellCounts.X << " x " << CellCounts.Y << " x " << CellCounts.Z
              << " (" << CellCounts.X * CellCounts.Y * CellCounts.Z << " cells)" << std::endl;

    int FilledCount = 0;
    for (int iz = MinCell.Z; iz <= MaxCell.Z; ++iz) {
        for (int iy = MinCell.Y; iy <= MaxCell.Y; ++iy) {
            for (int ix = MinCell.X; ix <= MaxCell.X; ++ix) {
                GS::Vector3i CellKey(ix, iy, iz);
                GS::AxisBox3d CellBounds = Grid.GetCellLocalBounds(CellKey);

                GS::AxisBox3d ExpandedBounds = CellBounds;
                GS::Vector3d BoundsSize = CellBounds.Max - CellBounds.Min;
                // ExpandedBounds.Expand(BoundsSize * 0.25);
                // if (!GS::TestBoxOverlap(ExpandedBounds, MeshInside))
                //     continue;

                // Compute BitGrid for this cell using sphere inside/outside
                auto TargetBG = SampleInsideBitGrid(CellBounds,
                    [&](const GS::Vector3d& P) { return MeshInside(P); });

                int TargetSetBits = TargetBG.CountSetBits();

                // skip mostly-empty cells
                if (TargetSetBits < SkipBitThreshold)
                    continue;

                // If all sample points are inside, fill with solid cell
                GS::ModelGridCell FillCell;
                if (TargetBG.CountSetBits() == BitGridNumBits) {
                    FillCell = GS::ModelGridCell::SolidCell();
                } else {
                    CellMatchResult Match = FindBestCellMatch_Precomputed(
                        TargetBG, CellTable);
                    FillCell = Match.Cell;
                    double FillPct = 100.0 * Match.MatchingBits / BitGridNumBits;
                    GS::ModelGridCellData_StandardRST::Parameters CellParams;
                    CellParams.Fields = FillCell.CellData;
                    // std::cout << "  Cell (" << ix << "," << iy << "," << iz
                    //           << ") type=" << (int)FillCell.CellType
                    //           << " dir=" << (int)CellParams.AxisDirection
                    //           << " rot=" << (int)CellParams.AxisRotation
                    //           << " match=" << FillPct << "%" << std::endl;
                }

                Editor.FillCell(CellKey, FillCell,
                    [](const GS::ModelGridCell&) { return true; },
                    [](const GS::ModelGridCell&, GS::ModelGridCell&) {});
                FilledCount++;
            }
        }
        std::cout << "  Z-layer " << (iz - MinCell.Z + 1) << "/" << CellCounts.Z << " done" << std::endl;
    }
    std::cout << "Filled " << FilledCount << " cells" << std::endl;

    // Generate the grid mesh
    GS::AxisBox3i OccupiedBounds = Grid.GetOccupiedRegionBounds(1);
    GS::AxisBox3d LocalBounds = Grid.GetCellLocalBounds(OccupiedBounds.Min);
    LocalBounds.Contain(Grid.GetCellLocalBounds(OccupiedBounds.Max));

    GS::DenseMeshBuilderFactory Factory;
    GS::ModelGridMeshCache MeshGen;
    MeshGen.Initialize(GS::Vector3d(CellSize, CellSize, CellSize), &Factory);
    MeshGen.UpdateInBounds(Grid, LocalBounds, [](GS::Vector2i) {});
    GS::DenseMeshCollector Collector;
    MeshGen.ExtractFullMesh(Collector);
    GS::DenseMesh GridMesh = Collector.AccumulatedMesh.ToDenseMesh();

    std::cout << "Grid mesh: " << GridMesh.GetVertexCount() << " vertices, "
              << GridMesh.GetTriangleCount() << " triangles" << std::endl;

    // Write meshes to OBJ
    GS::WriteMeshOBJ("../output/source_as_cells.obj", GridMesh, true);
    GS::WriteMeshOBJ("../output/source_mesh.obj", SourceMesh, true);

    return 0;
}
