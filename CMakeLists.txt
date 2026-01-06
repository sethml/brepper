cmake_minimum_required(VERSION 3.20)
project(brepper VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find dependencies
find_package(PCL 1.12 REQUIRED COMPONENTS common io filters segmentation surface)
find_package(OpenCASCADE REQUIRED)
find_package(Eigen3 REQUIRED)  # Remove version constraint for now
find_package(CLI11 CONFIG QUIET)

# OpenMP support (required for parallel processing)
if(APPLE)
    # Homebrew libomp paths
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
endif()
find_package(OpenMP REQUIRED)
message(STATUS "OpenMP found - parallel processing enabled")

# If CLI11 not found via package manager, fetch it
if(NOT CLI11_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        cli11
        GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
        GIT_TAG v2.4.2
    )
    FetchContent_MakeAvailable(cli11)
endif()

# Collect sources
file(GLOB_RECURSE SOURCES src/*.cpp)

# Main executable
add_executable(brepper ${SOURCES})

target_include_directories(brepper PRIVATE
    ${PCL_INCLUDE_DIRS}
    ${OpenCASCADE_INCLUDE_DIR}
    src
)

target_link_directories(brepper PRIVATE
    ${PCL_LIBRARY_DIRS}
)

# OpenCASCADE libraries needed
set(OCCT_LIBS
    TKernel TKMath TKG3d TKGeomBase TKBRep
    TKGeomAlgo TKTopAlgo TKShHealing TKFillet
    TKXSBase TKDESTEP
)

target_link_libraries(brepper PRIVATE
    ${PCL_LIBRARIES}
    ${OCCT_LIBS}
    Eigen3::Eigen
    CLI11::CLI11
    OpenMP::OpenMP_CXX
)

target_compile_definitions(brepper PRIVATE
    ${PCL_DEFINITIONS}
)

# Enable warnings
target_compile_options(brepper PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
)

# Install
install(TARGETS brepper RUNTIME DESTINATION bin)

# Optional: Tests
option(BUILD_TESTS "Build tests" OFF)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()