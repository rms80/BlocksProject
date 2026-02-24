#include <iostream>
#include <set>
#include <type_traits>

#include "Core/gs_parallel_api.h"
#include "ModelGrid/ModelGridMesher.h"
#include "ModelGrid/ModelGridCell.h"
#include "Math/GSTransformList.h"

#include "src/DenseMeshAPI.h"
#include "src/DummyParallelAPI.h"
#include "src/MeshCode.h"

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

template<typename CellType, int NumVertices>
void AccumulateDimensionCodes(GS::ModelGridMesher& Mesher,
    GS::DenseMeshBuilder& Builder, GS::ModelGridMesher::AppendCache& Cache,
    unsigned int TransX, unsigned int TransY, unsigned int TransZ,
    std::set<GS::MeshCode<NumVertices>>& Codes)
{
    constexpr unsigned int MaxDim = GS::ModelGridCellData_StandardRST::MaxDimension;
    constexpr unsigned int StartDim = MaxDim;
    GS::Vector3d CellDims(100.0, 100.0, 100.0);
    GS::AxisBox3d CellBounds(GS::Vector3d::Zero(), CellDims);
    GS::ModelGridMesher::CellMaterials Materials;

    constexpr unsigned int MaxAxis = GS::ModelGridCellData_StandardRST::MaxRotationAxis;
    constexpr unsigned int MaxRot = GS::ModelGridCellData_StandardRST::MaxRotationAngle;

    for (unsigned int Axis = 0; Axis <= MaxAxis; Axis++)
    {
        for (unsigned int Rot = 0; Rot <= MaxRot; Rot++)
        {
            for (unsigned int DimX = StartDim; DimX <= MaxDim; DimX++)
            {
                for (unsigned int DimY = StartDim; DimY <= MaxDim; DimY++)
                {
                    for (unsigned int DimZ = StartDim; DimZ <= MaxDim; DimZ++)
                    {
                        CellType Cell = CellType::GetDefaultCellParams();
                        Cell.Params.AxisDirection = Axis;
                        Cell.Params.AxisRotation = Rot;
                        Cell.Params.DimensionX = DimX;
                        Cell.Params.DimensionY = DimY;
                        Cell.Params.DimensionZ = DimZ;
                        Cell.Params.TranslateX = TransX;
                        Cell.Params.TranslateY = TransY;
                        Cell.Params.TranslateZ = TransZ;

                        GS::TransformListd Transforms;
                        GS::GetUnitCellTransform(Cell, CellDims, Transforms);

                        Builder.Reset();
                        Mesher.ResetAppendCache(Cache, false);
                        AppendCellMesh<CellType>(Mesher, CellBounds, Materials, Builder, Transforms, Cache);

                        GS::MeshCode<NumVertices> Code = GS::ComputeMeshCode<NumVertices>(Builder);
                        Codes.insert(Code);
                    }
                }
            }
        }
    }
}

template<typename CellType, int NumVertices>
void EnumerateCellCodes(GS::ModelGridMesher& Mesher,
    GS::DenseMeshBuilder& Builder, GS::ModelGridMesher::AppendCache& Cache,
    const char* CellName)
{
    std::set<GS::MeshCode<NumVertices>> UniqueCodes;

    AccumulateDimensionCodes<CellType, NumVertices>(Mesher, Builder, Cache, 0, 0, 0, UniqueCodes);
    std::cout << CellName << " unique codes (no translate): " << UniqueCodes.size() << std::endl;

    constexpr unsigned int MaxTrans = GS::ModelGridCellData_StandardRST::MaxTranslate;

    for (unsigned int TX = 0; TX <= MaxTrans; TX++)
    {
       for (unsigned int TY = 0; TY <= MaxTrans; TY++)
           AccumulateDimensionCodes<CellType, NumVertices>(Mesher, Builder, Cache, TX, TY, 0, UniqueCodes);
       std::cout << "  TranslateXY TX=" << TX << " codes=" << UniqueCodes.size() << std::endl;
    }
    std::cout << CellName << " after TranslateXY: " << UniqueCodes.size() << " unique codes" << std::endl;

    for (unsigned int TX = 0; TX <= MaxTrans; TX++)
    {
       for (unsigned int TZ = 0; TZ <= MaxTrans; TZ++)
           AccumulateDimensionCodes<CellType, NumVertices>(Mesher, Builder, Cache, TX, 0, TZ, UniqueCodes);
       std::cout << "  TranslateXZ TX=" << TX << " codes=" << UniqueCodes.size() << std::endl;
    }
    std::cout << CellName << " after TranslateXZ: " << UniqueCodes.size() << " unique codes" << std::endl;

    for (unsigned int TY = 0; TY <= MaxTrans; TY++)
    {
       for (unsigned int TZ = 0; TZ <= MaxTrans; TZ++)
           AccumulateDimensionCodes<CellType, NumVertices>(Mesher, Builder, Cache, 0, TY, TZ, UniqueCodes);
       std::cout << "  TranslateYZ TY=" << TY << " codes=" << UniqueCodes.size() << std::endl;
    }
    std::cout << CellName << " after TranslateYZ: " << UniqueCodes.size() << " unique codes" << std::endl;
}

template<typename CellType, int NumVertices>
void EnumerateFlipTest(GS::ModelGridMesher& Mesher,
    GS::DenseMeshBuilder& Builder, GS::ModelGridMesher::AppendCache& Cache,
    const char* CellName, unsigned int DimX, unsigned int DimY, unsigned int DimZ)
{
    GS::Vector3d CellDims(100.0, 100.0, 100.0);
    GS::AxisBox3d CellBounds(GS::Vector3d::Zero(), CellDims);
    GS::ModelGridMesher::CellMaterials Materials;

    constexpr unsigned int MaxAxis = GS::ModelGridCellData_StandardRST::MaxRotationAxis;
    constexpr unsigned int MaxRot = GS::ModelGridCellData_StandardRST::MaxRotationAngle;

    std::set<GS::MeshCode<NumVertices>> Codes;

    // Phase 1: enumerate all orientations/rotations with no flips
    for (unsigned int Axis = 0; Axis <= MaxAxis; Axis++)
    {
        for (unsigned int Rot = 0; Rot <= MaxRot; Rot++)
        {
            CellType Cell = CellType::GetDefaultCellParams();
            Cell.Params.AxisDirection = Axis;
            Cell.Params.AxisRotation = Rot;
            Cell.Params.DimensionX = DimX;
            Cell.Params.DimensionY = DimY;
            Cell.Params.DimensionZ = DimZ;

            GS::TransformListd Transforms;
            GS::GetUnitCellTransform(Cell, CellDims, Transforms);

            Builder.Reset();
            Mesher.ResetAppendCache(Cache, false);
            AppendCellMesh<CellType>(Mesher, CellBounds, Materials, Builder, Transforms, Cache);

            GS::MeshCode<NumVertices> Code = GS::ComputeMeshCode<NumVertices>(Builder);
            Codes.insert(Code);
        }
    }
    std::cout << CellName << " dims=(" << DimX << "," << DimY << "," << DimZ
              << ") orient/rot only (6x4=24): " << Codes.size() << " unique codes" << std::endl;

    // Phase 2: add flips (8 combinations of FlipX/Y/Z) for each orient/rot
    for (unsigned int FlipBits = 1; FlipBits < 8; FlipBits++)
    {
        unsigned int FX = (FlipBits >> 0) & 1;
        unsigned int FY = (FlipBits >> 1) & 1;
        unsigned int FZ = (FlipBits >> 2) & 1;

        for (unsigned int Axis = 0; Axis <= MaxAxis; Axis++)
        {
            for (unsigned int Rot = 0; Rot <= MaxRot; Rot++)
            {
                CellType Cell = CellType::GetDefaultCellParams();
                Cell.Params.AxisDirection = Axis;
                Cell.Params.AxisRotation = Rot;
                Cell.Params.DimensionX = DimX;
                Cell.Params.DimensionY = DimY;
                Cell.Params.DimensionZ = DimZ;
                Cell.Params.FlipX = FX;
                Cell.Params.FlipY = FY;
                Cell.Params.FlipZ = FZ;

                GS::TransformListd Transforms;
                GS::GetUnitCellTransform(Cell, CellDims, Transforms);

                Builder.Reset();
                Mesher.ResetAppendCache(Cache, false);
                AppendCellMesh<CellType>(Mesher, CellBounds, Materials, Builder, Transforms, Cache);

                GS::MeshCode<NumVertices> Code = GS::ComputeMeshCode<NumVertices>(Builder);
                Codes.insert(Code);
            }
        }
        std::cout << "  + flip(" << FX << "," << FY << "," << FZ << "): " << Codes.size() << " unique codes" << std::endl;
    }
    std::cout << CellName << " after all flips: " << Codes.size() << " unique codes" << std::endl;
}


int main()
{
    GS::Parallel::RegisterAPI(GSMakeUniquePtr<GS::DummyParallelAPI>());

    GS::ModelGridMesher Mesher;
    Mesher.Initialize(GS::Vector3d(100.0, 100.0, 100.0));

    GS::ModelGridMesher::AppendCache Cache;
    Mesher.InitAppendCache(Cache);
    GS::DenseMeshBuilder Builder;

    //EnumerateCellCodes<GS::MGCell_Ramp, 6>(Mesher, Builder, Cache, "Ramp");
    //EnumerateCellCodes<GS::MGCell_Slab, 8>(Mesher, Builder, Cache, "Slab");
    //EnumerateCellCodes<GS::MGCell_Corner, 4>(Mesher, Builder, Cache, "Corner");

    EnumerateFlipTest<GS::MGCell_CutCorner, 7>(Mesher, Builder, Cache, "CutCorner", 7, 7, 7);
    EnumerateFlipTest<GS::MGCell_CutCorner, 7>(Mesher, Builder, Cache, "CutCorner", 10, 7, 5);

    return 0;
}
