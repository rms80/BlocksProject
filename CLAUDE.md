# BlocksProject

## Build Scripts
- **run_current_demo.sh** — builds and runs the `cell_match_test` sample
- **viewer/build.sh** — builds the polyscope viewer. Use this instead of running cmake directly.

## Viewer
- Source: `viewer/main.cpp`, CMake: `viewer/CMakeLists.txt`
- Binary: `viewer/build/polyscope_viewer`
- Args: `--load-grid=path/to/grid.bin` to load a serialized ModelGrid
- Has cell selection via raycasting (ModelGridCollider), yellow wireframe box highlight
- Polyscope submodule at `viewer/polyscope/`

## Project Structure
- Submodules: `submodules/GradientspaceCore`, `GradientspaceIO`, `GradientspaceGrid`
- Shared src files: `src/` (DenseMeshAPI, DenseMeshAABBTree, GeometryUtils, BoxGenerator, SphereGenerator)
- Samples: `samples/` (enumerate_cells, visualize_cells, winding_test, cell_match_test, cut_cell)

## Build
- Platform: macOS, Unix Makefiles, Release config, C++20
- Main project build dir: `build/`
- Viewer build dir: `viewer/build/`
