cmake_minimum_required(VERSION 3.14)

project(polyscope_viewer)

set(CMAKE_CXX_STANDARD 20)

# Polyscope
add_subdirectory(polyscope)

# Gradientspace libraries (from parent project's submodules)
set(GS_SUBMODULES ${CMAKE_CURRENT_SOURCE_DIR}/../submodules)
add_subdirectory(${GS_SUBMODULES}/GradientspaceCore ${CMAKE_CURRENT_BINARY_DIR}/gscore_bin)
add_subdirectory(${GS_SUBMODULES}/GradientspaceIO ${CMAKE_CURRENT_BINARY_DIR}/gsio_bin)
add_subdirectory(${GS_SUBMODULES}/GradientspaceGrid ${CMAKE_CURRENT_BINARY_DIR}/gsgrid_bin)

# Shared src files from parent project
file(GLOB_RECURSE SRC_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../src/*.*")

add_executable(polyscope_viewer main.cpp ${SRC_FILES})

target_include_directories(polyscope_viewer PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../src
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(polyscope_viewer
    polyscope
    gradientspace_core
    gradientspace_io
    gradientspace_grid
)
