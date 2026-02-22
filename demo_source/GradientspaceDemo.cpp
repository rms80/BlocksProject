// GradientspaceDemo.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <filesystem>

#include "Core/TextIO.h"
#include "Core/BinaryIO.h"
#include "Mesh/DenseMesh.h"
#include "MeshIO/OBJReader.h"
#include "MeshIO/OBJWriter.h"
#include "MeshIO/STLReader.h"
#include "MeshIO/STLWriter.h"


#define ENABLE_GRADIENTSPACE_GRID
#ifdef ENABLE_GRADIENTSPACE_GRID
#include "ModelGrid/ModelGrid.h"
#include "ModelGrid/ModelGridCell.h"
#include "ModelGrid/ModelGridEditor.h"
#include "ModelGrid/ModelGridEditMachine.h"
#include "ModelGrid/ModelGridSerializer.h"
#include "ModelGrid/ModelGridMeshCache.h"
#endif

#include "src/DenseMeshAPI.h"
#include "src/DummyParallelAPI.h"

int main()
{
    std::cout << "Hello World!\n";

    GS::Parallel::RegisterAPI(GSMakeUniquePtr<GS::DummyParallelAPI>());

	auto WriteDenseMeshOBJ = [](GS::DenseMesh& Mesh, std::string filename)
	{
		GS::OBJFormatData WriteOBJData;
		GS::DenseMeshToOBJFormatData(Mesh, WriteOBJData);
		auto Writer = GS::FileTextWriter::OpenFile(filename);
		GS::OBJWriter::WriteOBJ(Writer, WriteOBJData);
		Writer.CloseFile();
	};

    GS::DenseMesh TmpMesh;

    std::filesystem::path outputDir = std::filesystem::current_path() / "output";
    std::filesystem::create_directories(outputDir);

    std::string testFilesPath = "test_files" + std::string(1, std::filesystem::path::preferred_separator);
    std::string writeFilesPath = outputDir.string() + std::string(1, std::filesystem::path::preferred_separator);

    // read OBJ file into a DenseMesh
    GS::OBJFormatData OBJData;
    bool bOBJReadOK = GS::OBJReader::ReadOBJ(testFilesPath+"bunny_open_200.obj", OBJData);
    std::cout << "OBJ Mesh Read ok: " << bOBJReadOK << std::endl;


#ifdef ENABLE_GRADIENTSPACE_GRID

    // create a ModelGrid object
    GS::ModelGrid Grid;
    Grid.Initialize(GS::Vector3d::One());


    // Create a ModelGridEditMachine, which is a staeful high-level "Editor" for
    // ModelGrid objects. It has a 'Cursor' which you can move around, while
    // doing operations. The code below fills a 2D rectangle.
    GS::ModelGridEditMachine EditMachine;
    EditMachine.Initialize(Grid);
    EditMachine.SetCurrentDrawCellType(GS::EModelGridCellType::Ramp_Parametric);
    EditMachine.SetActiveDrawPlaneNormal(GS::Vector3d::UnitZ());
    EditMachine.BeginSculptCells_Rect2D();
    EditMachine.SetInitialCellCursor(GS::Vector3i(-5, -5, 0), GS::Vector3d(-5,-5,0), GS::Vector3d::UnitZ());
    EditMachine.UpdateCellCursor(GS::Vector3i(5, 5, 0));    // fills 1x1 rectangle (-5 to +5, inclusive)
    EditMachine.EndCurrentInteraction();

    // ModelGridEditor is a class for doing low-level edits of the ModelGrid
    GS::ModelGridEditor Editor(Grid);

    int NumFilledCells = 0;
    // enumerate over all currently-filled cells
    Grid.EnumerateFilledCells([&](GS::Vector3i Key, const GS::ModelGridCell& CellInfo, GS::AxisBox3d LocalBounds) {
        NumFilledCells++;
        if ( abs(Key.X) < 3 && abs(Key.Y) < 3 )
        {
            // make a new "slab" cell type w/ default parameters
            GS::ModelGridCell SlabCell = GS::MakeDefaultCellFromType(GS::EModelGridCellType::Slab_Parametric);

            // use the Editor to fill the cell w/ the new type,
            // but 
            Editor.FillCell(Key, SlabCell,
                // this is a filter function that can be used to discard the edit 
                [](const GS::ModelGridCell& CurCell) { return true; },
                // this is a modifier function that can be used to transfer "current" cell params to the new cell
                [](const GS::ModelGridCell& CurCell, GS::ModelGridCell& NewCell) {
                    // do nothing for now - fully replace
                });
        }
        //Editor.PaintCell(GS::Vector3i::Zero(), GS::Color3b::Red());
    });

    // find the bounds of occupied cells in grid-space
    GS::AxisBox3i OccupiedBounds = Grid.GetOccupiedRegionBounds(1);
    // convert to local space
    GS::AxisBox3d LocalBounds = Grid.GetCellLocalBounds(OccupiedBounds.Min);
    LocalBounds.Contain(Grid.GetCellLocalBounds(OccupiedBounds.Max));

    GS::DenseMeshBuilderFactory Factory;
    GS::ModelGridMeshCache MeshGen;
    MeshGen.Initialize(GS::Vector3d::One(), &Factory);
    MeshGen.UpdateInBounds(Grid, LocalBounds, [](GS::Vector2i) {});
    GS::DenseMeshCollector Collector;
    MeshGen.ExtractFullMesh(Collector);   
    GS::DenseMesh CollectedMesh = Collector.AccumulatedMesh.ToDenseMesh();       
    WriteDenseMeshOBJ(CollectedMesh, writeFilesPath+"modelgrid_mesh.obj");



    // ModelGrid binary serialization
    // GS::MemorySerializer Serializer;
    // Serializer.BeginWrite();
    // bool bStoreOK = GS::ModelGridSerializer::Serialize(Grid, Serializer);
    // std::cout << "Grid write ok: " << bStoreOK << std::endl;
    // GS::ModelGrid RestoredGrid;
    // Serializer.BeginRead();
    // bool bRestoreOK = GS::ModelGridSerializer::Restore(RestoredGrid, Serializer);
    // std::cout << "Grid read ok: " << bStoreOK << std::endl;    


#endif
}

