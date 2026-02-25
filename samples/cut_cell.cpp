#include <iostream>

#include "Core/gs_parallel_api.h"

#include "src/DenseMeshAPI.h"
#include "src/DummyParallelAPI.h"
#include "src/GeometryUtils.h"
#include "src/BoxGenerator.h"


// Generate a cut-cell mesh by cutting the far corner off a cube [0,S]^3.
// Cut points are on the 3 edges from (S,S,S), at distance (step+1)*S/NumSteps from the far corner.
// The kept portion contains the origin (0,0,0).
// When a step is at max (NumSteps-1), the cut vertex coincides with a cube corner and is
// merged to avoid degenerate triangles. Face polygons are reduced accordingly.
GS::DenseMesh GenerateCutCornerParamMesh(double CellSize, int stepX, int stepY, int stepZ, int NumSteps)
{
    double dx = (stepX + 1) * CellSize / NumSteps;
    double dy = (stepY + 1) * CellSize / NumSteps;
    double dz = (stepZ + 1) * CellSize / NumSteps;
    double S = CellSize;

    GS::DenseMeshBuilder Builder;

    // Cube corner vertices (always present)
    int V3 = Builder.AppendVertex(GS::Vector3d(0, 0, 0));  // origin
    int V4 = Builder.AppendVertex(GS::Vector3d(S, 0, 0));  // +X
    int V5 = Builder.AppendVertex(GS::Vector3d(0, S, 0));  // +Y
    int V6 = Builder.AppendVertex(GS::Vector3d(0, 0, S));  // +Z
    int V7 = Builder.AppendVertex(GS::Vector3d(S, S, 0));  // XY
    int V8 = Builder.AppendVertex(GS::Vector3d(S, 0, S));  // XZ
    int V9 = Builder.AppendVertex(GS::Vector3d(0, S, S));  // YZ

    // Cut vertices: merge with coincident cube corner at max step
    int V0 = (stepX == NumSteps - 1) ? V9 : Builder.AppendVertex(GS::Vector3d(S-dx, S, S));
    int V1 = (stepY == NumSteps - 1) ? V8 : Builder.AppendVertex(GS::Vector3d(S, S-dy, S));
    int V2 = (stepZ == NumSteps - 1) ? V7 : Builder.AppendVertex(GS::Vector3d(S, S, S-dz));

    int GroupID = Builder.AllocateGroupID();

    // 7 face polygons (CW winding to match codebase convention)
    GS::FanTriangulate(Builder, GroupID, std::array{V0, V2, V1});               // Cut face
    GS::FanTriangulate(Builder, GroupID, std::array{V4, V8, V1, V2, V7});       // +X face
    GS::FanTriangulate(Builder, GroupID, std::array{V5, V7, V2, V0, V9});       // +Y face
    GS::FanTriangulate(Builder, GroupID, std::array{V6, V9, V0, V1, V8});       // +Z face
    GS::FanTriangulate(Builder, GroupID, std::array{V3, V5, V9, V6});           // -X face
    GS::FanTriangulate(Builder, GroupID, std::array{V3, V6, V8, V4});           // -Y face
    GS::FanTriangulate(Builder, GroupID, std::array{V3, V4, V7, V5});           // -Z face

    return Builder.ToDenseMesh();
}


// Generate a cut-cell mesh by cutting the far edge off a cube [0,S]^3.
// The far edge runs from (0,S,S) to (S,S,S) along X at max Y and max Z.
// stepTop controls how far down from y=S to cut on the "top" Y-edges (at z=S).
// stepFront controls how far back from z=S to cut on the "front" Z-edges (at y=S).
// The kept portion contains the origin (0,0,0).
GS::DenseMesh GenerateCutEdgeParamMesh(double CellSize, int stepTop, int stepFront, int NumSteps)
{
    double dt = (stepTop + 1) * CellSize / NumSteps;
    double df = (stepFront + 1) * CellSize / NumSteps;
    double S = CellSize;

    GS::DenseMeshBuilder Builder;

    // 6 cube corners (the 2 on the far edge are removed)
    int V0 = Builder.AppendVertex(GS::Vector3d(0, 0, 0));  // origin
    int V1 = Builder.AppendVertex(GS::Vector3d(S, 0, 0));  // +X
    int V2 = Builder.AppendVertex(GS::Vector3d(0, S, 0));  // +Y
    int V3 = Builder.AppendVertex(GS::Vector3d(0, 0, S));  // +Z
    int V4 = Builder.AppendVertex(GS::Vector3d(S, S, 0));  // XY
    int V5 = Builder.AppendVertex(GS::Vector3d(S, 0, S));  // XZ

    // 4 cut vertices on edges adjacent to the far edge; merge at max step
    int V6 = (stepTop == NumSteps - 1)   ? V3 : Builder.AppendVertex(GS::Vector3d(0, S-dt, S));  // left top
    int V7 = (stepTop == NumSteps - 1)   ? V5 : Builder.AppendVertex(GS::Vector3d(S, S-dt, S));  // right top
    int V8 = (stepFront == NumSteps - 1) ? V2 : Builder.AppendVertex(GS::Vector3d(0, S, S-df));  // left front
    int V9 = (stepFront == NumSteps - 1) ? V4 : Builder.AppendVertex(GS::Vector3d(S, S, S-df));  // right front

    int GroupID = Builder.AllocateGroupID();

    // 8 face polygons (CW winding to match codebase convention)
    GS::FanTriangulate(Builder, GroupID, std::array{V6, V8, V9, V7});       // Cut face
    GS::FanTriangulate(Builder, GroupID, std::array{V2, V4, V9, V8});       // +Y face
    GS::FanTriangulate(Builder, GroupID, std::array{V3, V6, V7, V5});       // +Z face
    GS::FanTriangulate(Builder, GroupID, std::array{V0, V2, V8, V6, V3});   // -X face
    GS::FanTriangulate(Builder, GroupID, std::array{V1, V5, V7, V9, V4});   // +X face
    GS::FanTriangulate(Builder, GroupID, std::array{V0, V3, V5, V1});       // -Y face
    GS::FanTriangulate(Builder, GroupID, std::array{V0, V1, V4, V2});       // -Z face

    return Builder.ToDenseMesh();
}


// Append a DenseMesh to a DenseMeshBuilder with a position offset
void AppendMeshWithOffset(GS::DenseMeshBuilder& Builder, const GS::DenseMesh& Mesh, const GS::Vector3d& Offset)
{
    int BaseVertex = (int)Builder.Vertices.size();
    for (int i = 0; i < Mesh.GetVertexCount(); ++i)
        Builder.AppendVertex(Mesh.GetPosition(i) + Offset);
    int GroupID = Builder.AllocateGroupID();
    for (int i = 0; i < Mesh.GetTriangleCount(); ++i) {
        GS::Index3i Tri = Mesh.GetTriangle(i);
        Builder.AppendTriangle(
            GS::Index3i(Tri.A + BaseVertex, Tri.B + BaseVertex, Tri.C + BaseVertex), GroupID);
    }
}


// Append axis indicator boxes at the given offset
void AppendAxisLines(GS::DenseMeshBuilder& AccumBuilder, double CellSize, const GS::Vector3d& Offset)
{
    double T = CellSize * 0.02;
    GS::DenseMeshBuilder AxisBuilder;

    GS::BoxGenerator XAxis;
    XAxis.Center = GS::Vector3d(CellSize * 0.5, 0, 0);
    XAxis.Dimensions = GS::Vector3d(CellSize, T, T);
    XAxis.Generate(AxisBuilder);

    GS::BoxGenerator YAxis;
    YAxis.Center = GS::Vector3d(0, CellSize * 0.5, 0);
    YAxis.Dimensions = GS::Vector3d(T, CellSize, T);
    YAxis.Generate(AxisBuilder);

    GS::BoxGenerator ZAxis;
    ZAxis.Center = GS::Vector3d(0, 0, CellSize * 0.5);
    ZAxis.Dimensions = GS::Vector3d(T, T, CellSize);
    ZAxis.Generate(AxisBuilder);

    AppendMeshWithOffset(AccumBuilder, AxisBuilder.ToDenseMesh(), Offset);
}


void visualize_cut_cells()
{
    double CellSize = 10.0;
    double Spacing = CellSize * 1.5;
    constexpr int NumSteps = 4;

    GS::DenseMeshBuilder AccumBuilder;
    int MeshCount = 0;

    for (int iz = 0; iz < NumSteps; ++iz) {
        for (int iy = 0; iy < NumSteps; ++iy) {
            for (int ix = 0; ix < NumSteps; ++ix) {

                GS::DenseMesh Mesh = GenerateCutCornerParamMesh(CellSize, ix, iy, iz, NumSteps);

                int NumDegen = GS::CountDegenerateTriangles(Mesh);
                if (NumDegen > 0)
                    std::cout << "  WARNING: cell (" << ix << "," << iy << "," << iz
                              << ") has " << NumDegen << " degenerate triangles" << std::endl;

                GS::Vector3d Offset(-ix * Spacing * 3.0, iy * Spacing, -iz * Spacing);
                AppendMeshWithOffset(AccumBuilder, Mesh, Offset);
                AppendAxisLines(AccumBuilder, CellSize, Offset);
                MeshCount++;
            }
        }
    }

    GS::DenseMesh Result = AccumBuilder.ToDenseMesh();
    std::cout << "visualize_cut_cells: " << MeshCount << " cell meshes, "
              << Result.GetVertexCount() << " vertices, "
              << Result.GetTriangleCount() << " triangles" << std::endl;

    GS::WriteMeshOBJ("../../output/cut_cell.obj", Result, true);
}


void visualize_edge_cut_cells()
{
    double CellSize = 10.0;
    double Spacing = CellSize * 1.5;
    constexpr int NumSteps = 4;

    GS::DenseMeshBuilder AccumBuilder;
    int MeshCount = 0;

    for (int ifront = 0; ifront < NumSteps; ++ifront) {
        for (int itop = 0; itop < NumSteps; ++itop) {

            GS::DenseMesh Mesh = GenerateCutEdgeParamMesh(CellSize, itop, ifront, NumSteps);

            int NumDegen = GS::CountDegenerateTriangles(Mesh);
            if (NumDegen > 0)
                std::cout << "  WARNING: edge cell (" << itop << "," << ifront
                          << ") has " << NumDegen << " degenerate triangles" << std::endl;

            GS::Vector3d Offset(-itop * Spacing * 3.0, 0, -ifront * Spacing);
            AppendMeshWithOffset(AccumBuilder, Mesh, Offset);
            AppendAxisLines(AccumBuilder, CellSize, Offset);
            MeshCount++;
        }
    }

    GS::DenseMesh Result = AccumBuilder.ToDenseMesh();
    std::cout << "visualize_edge_cut_cells: " << MeshCount << " cell meshes, "
              << Result.GetVertexCount() << " vertices, "
              << Result.GetTriangleCount() << " triangles" << std::endl;

    GS::WriteMeshOBJ("../../output/edge_cut_cell.obj", Result, true);
}


int main()
{
    GS::Parallel::RegisterAPI(GSMakeUniquePtr<GS::DummyParallelAPI>());

    visualize_cut_cells();
    visualize_edge_cut_cells();

    return 0;
}
