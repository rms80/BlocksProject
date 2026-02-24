#include <algorithm>
#include <cmath>
#include <iostream>

#include "Core/gs_parallel_api.h"
#include "ModelGrid/ModelGrid.h"
#include "ModelGrid/ModelGridEditor.h"
#include "ModelGrid/ModelGridMesher.h"
#include "ModelGrid/ModelGridMeshCache.h"
#include "ModelGrid/ModelGridCell.h"

#include "src/DenseMeshAPI.h"
#include "src/DummyParallelAPI.h"
#include "src/DenseMeshAABBTree.h"
#include "src/GeometryUtils.h"
#include "MeshIO/OBJReader.h"
#include "MeshIO/OBJFormatData.h"


int main()
{
    GS::Parallel::RegisterAPI(GSMakeUniquePtr<GS::DummyParallelAPI>());

    // Load source mesh
    GS::OBJFormatData OBJData;
    if (!GS::OBJReader::ReadOBJ("../input/stanford-bunny.obj", OBJData)) {
        std::cerr << "Failed to read input/stanford-bunny.obj" << std::endl;
        return 1;
    }
    GS::DenseMesh SourceMesh;
    GS::OBJFormatDataToDenseMesh(OBJData, SourceMesh);

    std::cout << "Source mesh: " << SourceMesh.GetVertexCount() << " vertices, "
              << SourceMesh.GetTriangleCount() << " triangles" << std::endl;

    // Build AABB tree for fast winding number
    GS::DenseMeshAABBTree SourceTree;
    SourceTree.Build(SourceMesh);

    // Compute cell size from mesh bounds
    constexpr int NumCellsPerDimension = 25;
    GS::AxisBox3d MeshBounds = GS::ComputeBounds(SourceMesh);
    double MaxDim = std::max({MeshBounds.DimensionX(), MeshBounds.DimensionY(), MeshBounds.DimensionZ()});
    double CellSize = MaxDim / NumCellsPerDimension;
    std::cout << "Cell size: " << CellSize << std::endl;

    GS::ModelGrid Grid;
    Grid.Initialize(GS::Vector3d(CellSize, CellSize, CellSize));
    GS::ModelGridEditor Editor(Grid);

    // Compute grid index range covering the mesh bounds
    bool bIsInGrid = false;
    GS::Vector3i MinCell = Grid.GetCellAtPosition(MeshBounds.Min, bIsInGrid);
    GS::Vector3i MaxCell = Grid.GetCellAtPosition(MeshBounds.Max, bIsInGrid);

    GS::Vector3i CellCounts = MaxCell - MinCell + GS::Vector3i(1, 1, 1);
    std::cout << "Cell range: " << CellCounts.X << " x " << CellCounts.Y << " x " << CellCounts.Z
              << " (" << CellCounts.X * CellCounts.Y * CellCounts.Z << " cells)" << std::endl;

    int FilledCount = 0;
    for (int iz = MinCell.Z; iz <= MaxCell.Z; ++iz) {
        for (int iy = MinCell.Y; iy <= MaxCell.Y; ++iy) {
            for (int ix = MinCell.X; ix <= MaxCell.X; ++ix) {
                GS::Vector3i CellKey(ix, iy, iz);
                GS::AxisBox3d CellBounds = Grid.GetCellLocalBounds(CellKey);
                GS::Vector3d Center = (CellBounds.Min + CellBounds.Max) * 0.5;

                if (SourceTree.FastWindingNumber(Center) > 0.5) {
                    Editor.FillCell(CellKey, GS::ModelGridCell::SolidCell(),
                        [](const GS::ModelGridCell&) { return true; },
                        [](const GS::ModelGridCell&, GS::ModelGridCell&) {});
                    FilledCount++;
                }
            }
        }
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

    GS::WriteMeshOBJ("../output/winding_test.obj", GridMesh, true);

    return 0;
}
