# Tests for brepper

# Use Catch2 for testing
include(FetchContent)
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.2
)
FetchContent_MakeAvailable(Catch2)

# Include Catch2 CMake extras for test discovery
list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

# Test data directory
set(TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/test_config.hpp.in"
    "${CMAKE_CURRENT_BINARY_DIR}/test_config.hpp"
)

# Collect all test source files
set(TEST_SOURCES
    test_stl_reader.cpp
    test_mesh_sampling.cpp
)

# Create test executable
add_executable(brepper_tests ${TEST_SOURCES})

target_include_directories(brepper_tests PRIVATE
    ${PCL_INCLUDE_DIRS}
    ${OpenCASCADE_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_directories(brepper_tests PRIVATE
    ${PCL_LIBRARY_DIRS}
)

# OpenCASCADE libraries needed (same as main executable)
set(OCCT_LIBS
    TKernel TKMath TKG3d TKGeomBase TKBRep
    TKGeomAlgo TKTopAlgo TKShHealing TKFillet
    TKXSBase TKDESTEP
)

# Link against brepper library sources (exclude main.cpp)
file(GLOB_RECURSE LIB_SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp)
list(FILTER LIB_SOURCES EXCLUDE REGEX "main\\.cpp$")

target_sources(brepper_tests PRIVATE ${LIB_SOURCES})

target_link_libraries(brepper_tests PRIVATE
    ${PCL_LIBRARIES}
    ${OCCT_LIBS}
    Eigen3::Eigen
    Catch2::Catch2WithMain
)

target_compile_definitions(brepper_tests PRIVATE
    ${PCL_DEFINITIONS}
)

# Register tests with CTest
catch_discover_tests(brepper_tests)
