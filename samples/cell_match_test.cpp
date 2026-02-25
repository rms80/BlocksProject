#include <array>
#include <cmath>
#include <fstream>
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
#include "ModelGrid/ModelGridCell_Extended.h"
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
#include "src/ConvexHull.h"
#include "Core/gs_serializer.h"
#include "ModelGrid/ModelGridSerializer.h"

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

template<>
GS::DenseMesh GenerateCellMesh<GS::MGCell_VariableCutCorner>(
    const GS::Vector3d& CellDims,
    const GS::MGCell_VariableCutCorner& Cell)
{
    GS::ModelGridMesher Mesher;
    Mesher.Initialize(CellDims);
    GS::DenseMeshBuilder Builder;

    GS::ModelGridCellData_StandardRST BaseParams;
    BaseParams.Params.Fields = Cell.Params.Fields;
    GS::TransformListd Transforms;
    GS::GetUnitCellTransform(BaseParams, CellDims, Transforms);

    GS::AxisBox3d CellBounds(GS::Vector3d::Zero(), CellDims);
    GS::ModelGridMesher::CellMaterials Materials;
    Mesher.AppendVariableCutCorner(CellBounds, Materials, Builder, Transforms,
        Cell.Params.ParamA, Cell.Params.ParamB, Cell.Params.ParamC);
    return Builder.ToDenseMesh();
}

template<>
GS::DenseMesh GenerateCellMesh<GS::MGCell_VariableCutEdge>(
    const GS::Vector3d& CellDims,
    const GS::MGCell_VariableCutEdge& Cell)
{
    GS::ModelGridMesher Mesher;
    Mesher.Initialize(CellDims);
    GS::DenseMeshBuilder Builder;

    GS::ModelGridCellData_StandardRST BaseParams;
    BaseParams.Params.Fields = Cell.Params.Fields;
    GS::TransformListd Transforms;
    GS::GetUnitCellTransform(BaseParams, CellDims, Transforms);

    GS::AxisBox3d CellBounds(GS::Vector3d::Zero(), CellDims);
    GS::ModelGridMesher::CellMaterials Materials;
    Mesher.AppendVariableCutEdge(CellBounds, Materials, Builder, Transforms,
        Cell.Params.ParamA, Cell.Params.ParamB);
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
    GS::ConvexHull<10> Hull;
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
    //const unsigned int TryDimensionVals[] = { 3, 5, 7, 9, 11, 13, 15 };
    //const unsigned int TryDimensionVals[] = { 1, 2, 3, 4, 5, 7, 9, 11, 13, 15 };
    const unsigned int TryDimensionVals[] = { 1, 3, 6, 9, 12, 15 };


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
                            Hull.BuildFromMesh(Mesh);
                            auto BG = SampleInsideBitGrid(LocalBounds,
                                [&](const GS::Vector3d& P) { return Hull.IsInside(P); });
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
    AddType(GS::MGCell_Corner{});
    //AddType(GS::MGCell_Pyramid{});
    //AddType(GS::MGCell_Peak{});
    //AddType(GS::MGCell_Cylinder{});
    //AddType(GS::MGCell_CutCorner{});

    //const unsigned int TryParamDimensionVals[] = { 3, 7, 11, 15 };
    const unsigned int TryParamDimensionVals[] = { 7, 15 };
    const unsigned int TryParamVals[] = { 1, 3, 6, 9, 12, 15 };
    constexpr size_t NumParamVals = std::size(TryParamVals);

    // AddVariableType: like AddType but also varies ParamA..ParamD from TryParamVals.
    // NumParams (1-4) controls how many of the extended params are varied;
    // the rest stay at their default value.
    auto AddVariableType = [&](auto CellTypeTag, int NumParams) {
        using CellType = decltype(CellTypeTag);
        const unsigned int* pRanges[4];
        size_t nRanges[4];
        CellType Defaults = CellType::GetDefaultCellParams();
        const unsigned int DefaultA[] = { (unsigned int)Defaults.Params.ParamA };
        const unsigned int DefaultB[] = { (unsigned int)Defaults.Params.ParamB };
        const unsigned int DefaultC[] = { (unsigned int)Defaults.Params.ParamC };
        const unsigned int DefaultD[] = { (unsigned int)Defaults.Params.ParamD };
        const unsigned int* DefaultVals[] = { DefaultA, DefaultB, DefaultC, DefaultD };
        for (int i = 0; i < 4; ++i) {
            if (i < NumParams) { pRanges[i] = TryParamVals; nRanges[i] = NumParamVals; }
            else               { pRanges[i] = DefaultVals[i]; nRanges[i] = 1; }
        }
        for (unsigned int Dir = 0; Dir < 6; ++Dir)
            for (unsigned int Rot = 0; Rot < 4; ++Rot)
                for (unsigned int DimX : TryParamDimensionVals)
                    for (unsigned int DimY : TryParamDimensionVals)
                        for (unsigned int DimZ : TryParamDimensionVals)
                            for (size_t iA = 0; iA < nRanges[0]; ++iA)
                                for (size_t iB = 0; iB < nRanges[1]; ++iB)
                                    for (size_t iC = 0; iC < nRanges[2]; ++iC)
                                        for (size_t iD = 0; iD < nRanges[3]; ++iD) {
                                            CellType TypedCell = CellType::GetDefaultCellParams();
                                            TypedCell.Params.AxisDirection = Dir;
                                            TypedCell.Params.AxisRotation = Rot;
                                            TypedCell.Params.DimensionX = DimX;
                                            TypedCell.Params.DimensionY = DimY;
                                            TypedCell.Params.DimensionZ = DimZ;
                                            TypedCell.Params.ParamA = pRanges[0][iA];
                                            TypedCell.Params.ParamB = pRanges[1][iB];
                                            TypedCell.Params.ParamC = pRanges[2][iC];
                                            TypedCell.Params.ParamD = pRanges[3][iD];
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

    AddVariableType(GS::MGCell_VariableCutCorner{}, 3);
    AddVariableType(GS::MGCell_VariableCutEdge{}, 2);

    std::cout << "  Filtered " << DuplicateCount << " duplicate BitGrids" << std::endl;
    return Table;
}


void FindBestCellMatches_Precomputed(
    const BitGrid& Target,
    const std::vector<PrecomputedCell>& Table,
    int N,
    std::vector<CellMatchResult>& Results)
{
    Results.clear();
    Results.reserve(N);

    int MinKept = 0;  // worst MatchingBits currently in Results

    for (const auto& Entry : Table) {
        int Matching = Entry.BG.CountMatchingBits(Target);
        if ((int)Results.size() < N || Matching > MinKept) {
            CellMatchResult R;
            R.Cell = Entry.Cell;
            R.MatchingBits = Matching;
            R.DifferentBits = BitGridNumBits - Matching;

            // Insert in sorted order (best first)
            auto it = Results.begin();
            while (it != Results.end() && it->MatchingBits >= Matching)
                ++it;
            Results.insert(it, R);

            if ((int)Results.size() > N)
                Results.resize(N);

            MinKept = Results.back().MatchingBits;

            if (Results[0].MatchingBits == BitGridNumBits && (int)Results.size() >= N)
                return;
        }
    }
}


// Local parameter-space refinement: for each of the N candidates, try nearby parameter
// values (DimensionX/Y/Z +/-1, extra params +/-1) and return the best overall match.
// Uses ConvexHull plane test instead of winding number for speed.
CellMatchResult RefineMatchLocally(
    const BitGrid& Target,
    const GS::Vector3d& CellDims,
    const std::vector<CellMatchResult>& Candidates)
{
    CellMatchResult Best;
    if (Candidates.empty())
        return Best;

    Best = Candidates[0];

    // If the best precomputed match is already >98%, just use it
    constexpr int HighMatchThreshold = (int)(BitGridNumBits * 0.98);
    if (Best.MatchingBits > HighMatchThreshold)
        return Best;

    // Early-out threshold: skip candidate if >15% worse than current best
    constexpr int EarlyOutGap = (int)(BitGridNumBits * 0.15);

    GS::AxisBox3d LocalBounds(GS::Vector3d::Zero(), CellDims);
    GS::ConvexHull<10> Hull;

    auto TryCell = [&](const GS::DenseMesh& Mesh) {
        Hull.BuildFromMesh(Mesh);
        auto BG = SampleInsideBitGrid(LocalBounds,
            [&](const GS::Vector3d& P) { return Hull.IsInside(P); });
        return BG.CountMatchingBits(Target);
    };

    auto ClampU = [](int v, int lo, int hi) -> unsigned int {
        return (unsigned int)std::max(lo, std::min(v, hi));
    };

    for (size_t ci = 0; ci < Candidates.size(); ++ci)
    {
        const auto& Candidate = Candidates[ci];
        GS::EModelGridCellType CellType = Candidate.Cell.CellType;

        // Skip empty/solid cells - no parameters to refine
        if (CellType == GS::EModelGridCellType::Empty || CellType == GS::EModelGridCellType::Filled)
            continue;

        // If this candidate is more than 15% worse than our current best, stop
        if (ci > 0 && Candidate.MatchingBits < Best.MatchingBits - EarlyOutGap)
            break;

        // Extract base parameters from the candidate
        GS::ModelGridCellData_StandardRST BaseRST;
        GS::InitializeSubCellFromGridCell(Candidate.Cell, BaseRST);
        int BaseDimX = (int)BaseRST.Params.DimensionX;
        int BaseDimY = (int)BaseRST.Params.DimensionY;
        int BaseDimZ = (int)BaseRST.Params.DimensionZ;

        // Determine extra params for variable types
        bool bIsVariableCutCorner = (CellType == GS::EModelGridCellType::VariableCutCorner_Parametric);
        bool bIsVariableCutEdge   = (CellType == GS::EModelGridCellType::VariableCutEdge_Parametric);
        bool bHasExtParams = bIsVariableCutCorner || bIsVariableCutEdge;

        int BaseParamA = 0, BaseParamB = 0, BaseParamC = 0;
        if (bHasExtParams) {
            GS::ModelGridCellData_StandardRST_Ext ExtParams;
            GS::InitializeSubCellFromGridCell(Candidate.Cell, ExtParams);
            BaseParamA = (int)ExtParams.Params.ParamA;
            BaseParamB = (int)ExtParams.Params.ParamB;
            BaseParamC = (int)ExtParams.Params.ParamC;
        }

        // Iterate over dimension range +/- 2
        constexpr int DimRange = 1;
        constexpr int ParamRange = 1;
        constexpr int MaxDim = 15;
        constexpr int MaxParam = 15;

        for (int ddx = -DimRange; ddx <= DimRange; ++ddx)
        for (int ddy = -DimRange; ddy <= DimRange; ++ddy)
        for (int ddz = -DimRange; ddz <= DimRange; ++ddz)
        {
            unsigned int DX = ClampU(BaseDimX + ddx, 0, MaxDim);
            unsigned int DY = ClampU(BaseDimY + ddy, 0, MaxDim);
            unsigned int DZ = ClampU(BaseDimZ + ddz, 0, MaxDim);

            if (bIsVariableCutCorner)
            {
                for (int dpA = -ParamRange; dpA <= ParamRange; ++dpA)
                for (int dpB = -ParamRange; dpB <= ParamRange; ++dpB)
                for (int dpC = -ParamRange; dpC <= ParamRange; ++dpC)
                {
                    GS::MGCell_VariableCutCorner Cell;
                    GS::InitializeSubCellFromGridCell(Candidate.Cell, Cell);
                    Cell.Params.DimensionX = DX;
                    Cell.Params.DimensionY = DY;
                    Cell.Params.DimensionZ = DZ;
                    Cell.Params.ParamA = ClampU(BaseParamA + dpA, 0, MaxParam);
                    Cell.Params.ParamB = ClampU(BaseParamB + dpB, 0, MaxParam);
                    Cell.Params.ParamC = ClampU(BaseParamC + dpC, 0, MaxParam);
                    GS::DenseMesh Mesh = GenerateCellMesh(CellDims, Cell);
                    int Matching = TryCell(Mesh);
                    if (Matching > Best.MatchingBits) {
                        GS::UpdateGridCellFromSubCell(Best.Cell, Cell);
                        Best.MatchingBits = Matching;
                        Best.DifferentBits = BitGridNumBits - Matching;
                    }
                }
            }
            else if (bIsVariableCutEdge)
            {
                for (int dpA = -ParamRange; dpA <= ParamRange; ++dpA)
                for (int dpB = -ParamRange; dpB <= ParamRange; ++dpB)
                {
                    GS::MGCell_VariableCutEdge Cell;
                    GS::InitializeSubCellFromGridCell(Candidate.Cell, Cell);
                    Cell.Params.DimensionX = DX;
                    Cell.Params.DimensionY = DY;
                    Cell.Params.DimensionZ = DZ;
                    Cell.Params.ParamA = ClampU(BaseParamA + dpA, 0, MaxParam);
                    Cell.Params.ParamB = ClampU(BaseParamB + dpB, 0, MaxParam);
                    GS::DenseMesh Mesh = GenerateCellMesh(CellDims, Cell);
                    int Matching = TryCell(Mesh);
                    if (Matching > Best.MatchingBits) {
                        GS::UpdateGridCellFromSubCell(Best.Cell, Cell);
                        Best.MatchingBits = Matching;
                        Best.DifferentBits = BitGridNumBits - Matching;
                    }
                }
            }
            else
            {
                // Standard types: just vary dimensions, dispatch by type
                auto TryStandardType = [&](auto CellTypeTag) {
                    using CT = decltype(CellTypeTag);
                    CT Cell;
                    GS::InitializeSubCellFromGridCell(Candidate.Cell, Cell);
                    Cell.Params.DimensionX = DX;
                    Cell.Params.DimensionY = DY;
                    Cell.Params.DimensionZ = DZ;
                    GS::DenseMesh Mesh = GenerateCellMesh<CT>(CellDims, Cell);
                    int Matching = TryCell(Mesh);
                    if (Matching > Best.MatchingBits) {
                        GS::UpdateGridCellFromSubCell(Best.Cell, Cell);
                        Best.MatchingBits = Matching;
                        Best.DifferentBits = BitGridNumBits - Matching;
                    }
                };

                switch (CellType) {
                    case GS::EModelGridCellType::Slab_Parametric:      TryStandardType(GS::MGCell_Slab{}); break;
                    case GS::EModelGridCellType::Ramp_Parametric:      TryStandardType(GS::MGCell_Ramp{}); break;
                    case GS::EModelGridCellType::Corner_Parametric:    TryStandardType(GS::MGCell_Corner{}); break;
                    case GS::EModelGridCellType::Pyramid_Parametric:   TryStandardType(GS::MGCell_Pyramid{}); break;
                    case GS::EModelGridCellType::Peak_Parametric:      TryStandardType(GS::MGCell_Peak{}); break;
                    case GS::EModelGridCellType::Cylinder_Parametric:  TryStandardType(GS::MGCell_Cylinder{}); break;
                    case GS::EModelGridCellType::CutCorner_Parametric: TryStandardType(GS::MGCell_CutCorner{}); break;
                    default: break;
                }
            }
        }
    }

    return Best;
}


int main()
{
    GS::Parallel::RegisterAPI(GSMakeUniquePtr<GS::DummyParallelAPI>());

    GS::DenseMesh SourceMesh;

    std::string InputFilePath = "../input/stanford-bunny.obj";
    //std::string InputFilePath = "../input/Bauhaus_main.obj";

    // Load source mesh
    GS::OBJFormatData OBJData;
    if (!GS::OBJReader::ReadOBJ(InputFilePath, OBJData)) {
        std::cerr << "Failed to read " << InputFilePath << std::endl;
        return 1;
    }
    GS::OBJFormatDataToDenseMesh(OBJData, SourceMesh);

    // GS::DenseMeshBuilder Builder;
    // GS::BoxGenerator BoxGen;
    // BoxGen.Center = GS::Vector3d(2.0, 5.0, 7.0);
    // BoxGen.Dimensions = GS::Vector3d(100.0, 100.0, 100.0);
    // BoxGen.Generate(Builder);
    // SourceMesh = Builder.ToDenseMesh();

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
    std::vector<CellMatchResult> CellMatches;

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
                // if (TargetSetBits < SkipBitThreshold)
                //     continue;

                // If all sample points are inside, fill with solid cell
                GS::ModelGridCell FillCell;
                if (TargetBG.CountSetBits() == BitGridNumBits) {
                    FillCell = GS::ModelGridCell::SolidCell();
                } else {
                    FindBestCellMatches_Precomputed(
                        TargetBG, CellTable, 3, CellMatches);
                    CellMatchResult Refined = RefineMatchLocally(
                        TargetBG, CellDims, CellMatches);
                    FillCell = Refined.Cell;
                    double FillPct = 100.0 * Refined.MatchingBits / BitGridNumBits;
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

    // Write binary grid file
    {
        GS::MemorySerializer Serializer;
        Serializer.BeginWrite();
        bool bStoreOK = GS::ModelGridSerializer::Serialize(Grid, Serializer);
        if (bStoreOK)
        {
            size_t NumBytes = 0;
            const uint8_t* Buffer = Serializer.GetBuffer(NumBytes);
            std::ofstream outFile("../output/source_as_cells.grid", std::ios::binary);
            outFile.write(reinterpret_cast<const char*>(Buffer), NumBytes);
            outFile.close();
            std::cout << "Wrote ../output/source_as_cells.grid (" << NumBytes << " bytes)" << std::endl;
        }
        else
        {
            std::cout << "Grid serialization failed!" << std::endl;
        }
    }

    return 0;
}
