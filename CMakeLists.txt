cmake_minimum_required(VERSION 3.22)

project(gradientspace_demo)

# this gets rid of MinSizeRel and RelWithDebInfo configurations...even though they are still listed in CMake GUI...
#set(CMAKE_CONFIGURATION_TYPES "Debug;Release")

# this sets CMAKE_INSTALL_PREFIX automatically, otherwise it defaults to ProgramFiles(x86)...
#IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
#  SET(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/install CACHE PATH "install folder is not used?" FORCE)
#ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

#use C++20
set (CMAKE_CXX_STANDARD 20)


file(GLOB_RECURSE DEMO_SOURCE_FILES DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/demo_source/*.*")
file(GLOB_RECURSE SRC_FILES DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/*.*")
set(MODULE_FILES ${DEMO_SOURCE_FILES} ${SRC_FILES})

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/submodules/GradientspaceCore ${CMAKE_CURRENT_BINARY_DIR}/gscore_bin )
add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/submodules/GradientspaceIO ${CMAKE_CURRENT_BINARY_DIR}/gsio_bin )
add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/submodules/GradientspaceGrid ${CMAKE_CURRENT_BINARY_DIR}/gsgrid_bin )

add_executable(gradientspace_demo ${MODULE_FILES})

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT gradientspace_demo)

target_include_directories(gradientspace_demo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(gradientspace_demo gradientspace_core)
target_link_libraries(gradientspace_demo gradientspace_io)
target_link_libraries(gradientspace_demo gradientspace_grid)


# enumerate_cells sample
add_executable(enumerate_cells samples/enumerate_cells.cpp ${SRC_FILES})
target_include_directories(enumerate_cells PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(enumerate_cells gradientspace_core)
target_link_libraries(enumerate_cells gradientspace_io)
target_link_libraries(enumerate_cells gradientspace_grid)


# visualize_cells sample
add_executable(visualize_cells samples/visualize_cells.cpp ${SRC_FILES})
target_include_directories(visualize_cells PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(visualize_cells gradientspace_core)
target_link_libraries(visualize_cells gradientspace_io)
target_link_libraries(visualize_cells gradientspace_grid)


# winding_test sample
add_executable(winding_test samples/winding_test.cpp ${SRC_FILES})
target_include_directories(winding_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(winding_test gradientspace_core)
target_link_libraries(winding_test gradientspace_io)
target_link_libraries(winding_test gradientspace_grid)


# cell_match_test sample
add_executable(cell_match_test samples/cell_match_test.cpp ${SRC_FILES})
target_include_directories(cell_match_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(cell_match_test gradientspace_core)
target_link_libraries(cell_match_test gradientspace_io)
target_link_libraries(cell_match_test gradientspace_grid)


# cut_cell sample
add_executable(cut_cell samples/cut_cell.cpp ${SRC_FILES})
target_include_directories(cut_cell PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(cut_cell gradientspace_core)
target_link_libraries(cut_cell gradientspace_io)
target_link_libraries(cut_cell gradientspace_grid)


# from https://stackoverflow.com/questions/14089284/copy-all-dlls-that-an-executable-links-to-to-the-executable-directory
# this copies all the dependent DLLs to the folder the exe was generaed in
# however it doesn't copy PDBs...

# The following variable is defined only on DLL systems
if (CMAKE_IMPORT_LIBRARY_SUFFIX)
  add_custom_command(
    TARGET gradientspace_demo POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_RUNTIME_DLLS:gradientspace_demo> $<TARGET_FILE_DIR:gradientspace_demo>
    COMMAND_EXPAND_LISTS
  )
endif ()

# create a symlink to the executable in the source directory
add_custom_command(
  TARGET gradientspace_demo POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E rm -f ${CMAKE_CURRENT_SOURCE_DIR}/gradientspace_demo
  COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE:gradientspace_demo> ${CMAKE_CURRENT_SOURCE_DIR}/gradientspace_demo
)

