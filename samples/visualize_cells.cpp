#include <iostream>
#include <type_traits>

#include "Core/gs_parallel_api.h"
#include "ModelGrid/ModelGridMesher.h"
#include "ModelGrid/ModelGridCell.h"
#include "Math/GSTransformList.h"

#include "src/DenseMeshAPI.h"
#include "src/DummyParallelAPI.h"
#include "src/GeometryUtils.h"
#include "src/BoxGenerator.h"


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


void visualize_dim()
{
    double CellSize = 10.0;
    GS::Vector3d CellDims(CellSize, CellSize, CellSize);
    double Spacing = CellSize * 1.5;

    const unsigned int DimVals[] = { 0, 5, 10, 15 };
    constexpr int NumSteps = 4;

    GS::DenseMeshBuilder AccumBuilder;
    int MeshCount = 0;

    for (int iz = 0; iz < NumSteps; ++iz) {
        for (int iy = 0; iy < NumSteps; ++iy) {
            for (int ix = 0; ix < NumSteps; ++ix) {
                GS::MGCell_Corner Cell = GS::MGCell_Corner::GetDefaultCellParams();
                Cell.Params.DimensionX = DimVals[ix];
                Cell.Params.DimensionY = DimVals[iy];
                Cell.Params.DimensionZ = DimVals[iz];

                GS::DenseMesh Mesh = GenerateCellMesh<GS::MGCell_Corner>(CellDims, Cell);

                GS::Vector3d Offset(-ix * Spacing * 3.0, iy * Spacing, -iz * Spacing);
                AppendMeshWithOffset(AccumBuilder, Mesh, Offset);
                AppendAxisLines(AccumBuilder, CellSize, Offset);
                MeshCount++;
            }
        }
    }

    GS::DenseMesh Result = AccumBuilder.ToDenseMesh();
    std::cout << "visualize_dim: " << MeshCount << " cell meshes, "
              << Result.GetVertexCount() << " vertices, "
              << Result.GetTriangleCount() << " triangles" << std::endl;

    GS::WriteMeshOBJ("../output/visualize_dim.obj", Result, true);
}


void visualize_orient()
{
    double CellSize = 10.0;
    GS::Vector3d CellDims(CellSize, CellSize, CellSize);
    double Spacing = CellSize * 1.5;

    GS::DenseMeshBuilder AccumBuilder;
    int MeshCount = 0;

    for (unsigned int Rot = 0; Rot < 4; ++Rot) {
        for (unsigned int Dir = 0; Dir < 6; ++Dir) {
            GS::MGCell_Corner Cell = GS::MGCell_Corner::GetDefaultCellParams();
            Cell.Params.DimensionX = 5;
            Cell.Params.DimensionY = 5;
            Cell.Params.DimensionZ = 9;
            Cell.Params.AxisDirection = Dir;
            Cell.Params.AxisRotation = Rot;

            GS::DenseMesh Mesh = GenerateCellMesh<GS::MGCell_Corner>(CellDims, Cell);

            GS::Vector3d Offset(Dir * Spacing, -((double)Rot) * Spacing, 0);
            AppendMeshWithOffset(AccumBuilder, Mesh, Offset);
            AppendAxisLines(AccumBuilder, CellSize, Offset);
            MeshCount++;
        }
    }

    GS::DenseMesh Result = AccumBuilder.ToDenseMesh();
    std::cout << "visualize_orient: " << MeshCount << " cell meshes, "
              << Result.GetVertexCount() << " vertices, "
              << Result.GetTriangleCount() << " triangles" << std::endl;

    GS::WriteMeshOBJ("../output/visualize_orient.obj", Result, true);
}


int main()
{
    GS::Parallel::RegisterAPI(GSMakeUniquePtr<GS::DummyParallelAPI>());

    visualize_dim();
    visualize_orient();

    return 0;
}
