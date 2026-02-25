#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <algorithm>

#include "polyscope/polyscope.h"
#include "polyscope/options.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/view.h"
#include "imgui.h"

#include "Core/gs_parallel_api.h"
#include "Core/gs_serializer.h"
#include "Mesh/DenseMesh.h"
#include "Math/GSRay3.h"
#include "ModelGrid/ModelGrid.h"
#include "ModelGrid/ModelGridCell.h"
#include "ModelGrid/ModelGridEditor.h"
#include "ModelGrid/ModelGridEditMachine.h"
#include "ModelGrid/ModelGridMeshCache.h"
#include "ModelGrid/ModelGridSerializer.h"
#include "ModelGrid/ModelGridCollision.h"

#include "src/DenseMeshAPI.h"
#include "src/DummyParallelAPI.h"


// --- Global viewer state ---
struct ViewerState
{
    GS::ModelGrid* Grid = nullptr;
    GS::ModelGridCollider Collider;

    bool HasSelection = false;
    GS::Vector3i SelectedCellKey;

    bool IsolateMode = false;
    bool SlicePlaneEnabled = false;
    int SlicePlaneAxis = 2;       // 0=X, 1=Y, 2=Z
    float SlicePlaneOffset = 0.0f;
};
static ViewerState gState;


polyscope::SurfaceMesh* RegisterDenseMesh(const std::string& Name, const GS::DenseMesh& Mesh)
{
    int NumVerts = Mesh.GetVertexCount();
    int NumTris = Mesh.GetTriangleCount();

    std::vector<std::array<double, 3>> Vertices(NumVerts);
    for (int i = 0; i < NumVerts; i++) {
        const GS::Vector3d& P = Mesh.GetPosition(i);
        Vertices[i] = {P.X, P.Y, P.Z};
    }

    std::vector<std::array<uint32_t, 3>> Faces(NumTris);
    for (int i = 0; i < NumTris; i++) {
        const GS::Index3i& T = Mesh.GetTriangle(i);
        Faces[i] = {(uint32_t)T.A, (uint32_t)T.B, (uint32_t)T.C};
    }

    return polyscope::registerSurfaceMesh(Name, Vertices, Faces);
}


void UpdateSelectionViz()
{
    const std::string BoxName = "selection_box";

    if (!gState.HasSelection) {
        polyscope::removeCurveNetwork(BoxName, false);
        return;
    }

    GS::AxisBox3d Bounds = gState.Grid->GetCellLocalBounds(gState.SelectedCellKey);
    GS::Vector3d Lo = Bounds.Min;
    GS::Vector3d Hi = Bounds.Max;

    // 8 corners of the box
    std::vector<std::array<double, 3>> Nodes = {
        {Lo.X, Lo.Y, Lo.Z},  // 0
        {Hi.X, Lo.Y, Lo.Z},  // 1
        {Hi.X, Hi.Y, Lo.Z},  // 2
        {Lo.X, Hi.Y, Lo.Z},  // 3
        {Lo.X, Lo.Y, Hi.Z},  // 4
        {Hi.X, Lo.Y, Hi.Z},  // 5
        {Hi.X, Hi.Y, Hi.Z},  // 6
        {Lo.X, Hi.Y, Hi.Z},  // 7
    };

    // 12 edges of the box
    std::vector<std::array<size_t, 2>> Edges = {
        {0,1}, {1,2}, {2,3}, {3,0},   // bottom face
        {4,5}, {5,6}, {6,7}, {7,4},   // top face
        {0,4}, {1,5}, {2,6}, {3,7},   // vertical edges
    };

    auto* CN = polyscope::registerCurveNetwork(BoxName, Nodes, Edges);
    CN->setColor({1.0f, 1.0f, 0.0f});
    CN->setRadius(0.002f, true);
}


void RebuildViewMesh();   // forward declaration

void ViewerCallback()
{
    static constexpr float DragThreshold = 3.0f;
    static bool MouseWasPressed = false;
    static glm::vec2 MouseDownPos{0, 0};

    ImGuiIO& IO = ImGui::GetIO();

    // Track mouse-down position
    if (!IO.WantCaptureMouse && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
        MouseWasPressed = true;
        MouseDownPos = {IO.MousePos.x, IO.MousePos.y};
    }

    // On release, only select if the mouse didn't drag
    if (MouseWasPressed && ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
        MouseWasPressed = false;
        glm::vec2 ReleasePos{IO.MousePos.x, IO.MousePos.y};
        float Dist = glm::length(ReleasePos - MouseDownPos);

        if (Dist < DragThreshold) {
            glm::vec3 RayDir = polyscope::view::screenCoordsToWorldRay(ReleasePos);
            glm::vec3 CamPos = polyscope::view::getCameraWorldPosition();

            GS::Vector3d Origin(CamPos.x, CamPos.y, CamPos.z);
            GS::Vector3d Direction(RayDir.x, RayDir.y, RayDir.z);
            Direction.Normalize();
            GS::Ray3d Ray(Origin, Direction);

            double HitParam = 0;
            GS::Vector3d HitNormal;
            GS::Vector3i HitCellKey;

            if (gState.Collider.FindNearestHitCell(Ray, HitParam, HitNormal, HitCellKey)) {
                gState.HasSelection = true;
                gState.SelectedCellKey = HitCellKey;
                std::cout << "Selected cell: (" << HitCellKey.X << ", "
                          << HitCellKey.Y << ", " << HitCellKey.Z << ")" << std::endl;
            } else {
                gState.HasSelection = false;
            }

            UpdateSelectionViz();

            if (gState.IsolateMode) {
                RebuildViewMesh();
            }
        }
    }

    // Show selection info in a small UI panel
    ImGui::Begin("Selection");
    if (gState.HasSelection) {
        ImGui::Text("Cell: (%d, %d, %d)",
            gState.SelectedCellKey.X, gState.SelectedCellKey.Y, gState.SelectedCellKey.Z);

        if (ImGui::Button(gState.IsolateMode ? "Show All" : "Isolate Cell")) {
            gState.IsolateMode = !gState.IsolateMode;
            RebuildViewMesh();
        }
    } else {
        ImGui::Text("No cell selected");
        if (gState.IsolateMode) {
            gState.IsolateMode = false;
            RebuildViewMesh();
        }
    }

    ImGui::Separator();
    ImGui::Text("Slice Plane");

    bool SliceChanged = false;
    SliceChanged |= ImGui::Checkbox("Enable Slice", &gState.SlicePlaneEnabled);

    if (gState.SlicePlaneEnabled) {
        SliceChanged |= ImGui::RadioButton("X", &gState.SlicePlaneAxis, 0); ImGui::SameLine();
        SliceChanged |= ImGui::RadioButton("Y", &gState.SlicePlaneAxis, 1); ImGui::SameLine();
        SliceChanged |= ImGui::RadioButton("Z", &gState.SlicePlaneAxis, 2);

        GS::AxisBox3i OccBounds = gState.Grid->GetOccupiedRegionBounds(1);
        GS::AxisBox3d LocalMin = gState.Grid->GetCellLocalBounds(OccBounds.Min);
        GS::AxisBox3d LocalMax = gState.Grid->GetCellLocalBounds(OccBounds.Max);
        float SliderMin = 0.0f, SliderMax = 0.0f;
        if (gState.SlicePlaneAxis == 0) {
            SliderMin = (float)LocalMin.Min.X; SliderMax = (float)LocalMax.Max.X;
        } else if (gState.SlicePlaneAxis == 1) {
            SliderMin = (float)LocalMin.Min.Y; SliderMax = (float)LocalMax.Max.Y;
        } else {
            SliderMin = (float)LocalMin.Min.Z; SliderMax = (float)LocalMax.Max.Z;
        }
        gState.SlicePlaneOffset = std::clamp(gState.SlicePlaneOffset, SliderMin, SliderMax);
        SliceChanged |= ImGui::SliderFloat("Offset", &gState.SlicePlaneOffset, SliderMin, SliderMax);
    }

    if (SliceChanged && !gState.IsolateMode) {
        RebuildViewMesh();
    }

    ImGui::End();
}


GS::DenseMesh MeshModelGrid(const GS::ModelGrid& Grid)
{
    GS::AxisBox3i OccupiedBounds = Grid.GetOccupiedRegionBounds(1);
    GS::AxisBox3d LocalBounds = Grid.GetCellLocalBounds(OccupiedBounds.Min);
    LocalBounds.Contain(Grid.GetCellLocalBounds(OccupiedBounds.Max));

    GS::DenseMeshBuilderFactory Factory;
    GS::ModelGridMeshCache MeshGen;
    MeshGen.Initialize(Grid.GetCellDimensions(), &Factory);
    MeshGen.UpdateInBounds(Grid, LocalBounds, [](GS::Vector2i) {});
    GS::DenseMeshCollector Collector;
    MeshGen.ExtractFullMesh(Collector);
    return Collector.AccumulatedMesh.ToDenseMesh();
}


GS::DenseMesh MeshSingleCell(const GS::ModelGrid& Grid, GS::Vector3i CellKey)
{
    bool bIsInGrid = false;
    GS::ModelGridCell CellInfo = Grid.GetCellInfo(CellKey, bIsInGrid);
    if (!bIsInGrid) {
        return GS::DenseMesh();
    }

    GS::ModelGrid TempGrid;
    TempGrid.Initialize(Grid.GetCellDimensions());

    GS::ModelGridEditor Editor(TempGrid);
    Editor.FillCell(CellKey, CellInfo,
        [](const GS::ModelGridCell&) { return true; },
        [](const GS::ModelGridCell&, GS::ModelGridCell&) {});

    return MeshModelGrid(TempGrid);
}


GS::DenseMesh MeshSlicedGrid(const GS::ModelGrid& Grid, int Axis, float Offset)
{
    GS::ModelGrid TempGrid;
    TempGrid.Initialize(Grid.GetCellDimensions());
    GS::ModelGridEditor Editor(TempGrid);

    // Need a non-const reference for EnumerateFilledCells with bounds
    GS::ModelGrid& MutableGrid = const_cast<GS::ModelGrid&>(Grid);
    MutableGrid.EnumerateFilledCells([&](GS::Vector3i Key, const GS::ModelGridCell& CellInfo, GS::AxisBox3d LocalBounds) {
        GS::Vector3d Center = (LocalBounds.Min + LocalBounds.Max) * 0.5;
        double CenterVal = (Axis == 0) ? Center.X : (Axis == 1) ? Center.Y : Center.Z;
        if (CenterVal <= (double)Offset) {
            Editor.FillCell(Key, CellInfo,
                [](const GS::ModelGridCell&) { return true; },
                [](const GS::ModelGridCell&, GS::ModelGridCell&) {});
        }
    });

    return MeshModelGrid(TempGrid);
}


void RebuildViewMesh()
{
    GS::DenseMesh Mesh;

    if (gState.IsolateMode && gState.HasSelection) {
        Mesh = MeshSingleCell(*gState.Grid, gState.SelectedCellKey);
    } else if (gState.SlicePlaneEnabled) {
        Mesh = MeshSlicedGrid(*gState.Grid, gState.SlicePlaneAxis, gState.SlicePlaneOffset);
    } else {
        Mesh = MeshModelGrid(*gState.Grid);
    }

    RegisterDenseMesh("ModelGrid", Mesh)->setEdgeWidth(1.0);
}


bool LoadModelGridFromFile(const std::string& FilePath, GS::ModelGrid& Grid)
{
    std::ifstream File(FilePath, std::ios::binary | std::ios::ate);
    if (!File.is_open()) {
        std::cerr << "Failed to open file: " << FilePath << std::endl;
        return false;
    }

    size_t FileSize = File.tellg();
    File.seekg(0, std::ios::beg);

    GS::MemorySerializer Serializer;
    Serializer.InitializeMemory(FileSize);
    size_t BufferSize = 0;
    uint8_t* Buffer = Serializer.GetWritableBuffer(BufferSize);
    if (!File.read(reinterpret_cast<char*>(Buffer), FileSize)) {
        std::cerr << "Failed to read file: " << FilePath << std::endl;
        return false;
    }
    File.close();

    Serializer.BeginRead();
    if (!GS::ModelGridSerializer::Restore(Grid, Serializer)) {
        std::cerr << "Failed to deserialize ModelGrid from: " << FilePath << std::endl;
        return false;
    }

    return true;
}


GS::ModelGrid CreateSampleGrid()
{
    GS::ModelGrid Grid;
    Grid.Initialize(GS::Vector3d::One());

    GS::ModelGridEditMachine EditMachine;
    EditMachine.Initialize(Grid);
    EditMachine.SetCurrentDrawCellType(GS::EModelGridCellType::Ramp_Parametric);
    EditMachine.SetActiveDrawPlaneNormal(GS::Vector3d::UnitZ());
    EditMachine.BeginSculptCells_Rect2D();
    EditMachine.SetInitialCellCursor(GS::Vector3i(-5, -5, 0), GS::Vector3d(-5, -5, 0), GS::Vector3d::UnitZ());
    EditMachine.UpdateCellCursor(GS::Vector3i(5, 5, 0));
    EditMachine.EndCurrentInteraction();

    GS::ModelGridEditor Editor(Grid);
    Grid.EnumerateFilledCells([&](GS::Vector3i Key, const GS::ModelGridCell& CellInfo, GS::AxisBox3d LocalBounds) {
        if (abs(Key.X) < 3 && abs(Key.Y) < 3) {
            GS::ModelGridCell SlabCell = GS::MakeDefaultCellFromType(GS::EModelGridCellType::Slab_Parametric);
            Editor.FillCell(Key, SlabCell,
                [](const GS::ModelGridCell&) { return true; },
                [](const GS::ModelGridCell&, GS::ModelGridCell&) {});
        }
    });

    return Grid;
}


int main(int argc, char* argv[])
{
    GS::Parallel::RegisterAPI(GSMakeUniquePtr<GS::DummyParallelAPI>());

    GS::ModelGrid Grid;

    std::string GridFilePath;
    const std::string LoadGridPrefix = "--load-grid=";
    for (int i = 1; i < argc; i++) {
        std::string Arg = argv[i];
        if (Arg.rfind(LoadGridPrefix, 0) == 0) {
            GridFilePath = Arg.substr(LoadGridPrefix.size());
        }
    }

    if (!GridFilePath.empty()) {
        std::cout << "Loading ModelGrid from: " << GridFilePath << std::endl;
        if (!LoadModelGridFromFile(GridFilePath, Grid)) {
            return 1;
        }
        std::cout << "ModelGrid loaded successfully." << std::endl;
    } else {
        std::cout << "No --load-grid specified, creating sample grid." << std::endl;
        Grid = CreateSampleGrid();
    }

    // Build collision data
    gState.Grid = &Grid;
    gState.Collider.Initialize(Grid);
    GS::AxisBox3i OccupiedBounds = Grid.GetOccupiedRegionBounds(1);
    GS::AxisBox3d CollisionBounds = Grid.GetCellLocalBounds(OccupiedBounds.Min);
    CollisionBounds.Contain(Grid.GetCellLocalBounds(OccupiedBounds.Max));
    gState.Collider.UpdateInBounds(Grid, CollisionBounds);

    GS::DenseMesh GridMesh = MeshModelGrid(Grid);
    std::cout << "Grid mesh: " << GridMesh.GetVertexCount() << " vertices, "
              << GridMesh.GetTriangleCount() << " triangles" << std::endl;

    polyscope::options::uiScale = 2.0f / 3.0f;
    polyscope::init();
    RegisterDenseMesh("ModelGrid", GridMesh)->setEdgeWidth(1.0);
    polyscope::state::userCallback = ViewerCallback;
    polyscope::show();

    return 0;
}
