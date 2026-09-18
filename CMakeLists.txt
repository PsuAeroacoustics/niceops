cmake_minimum_required(VERSION 3.16)

project(NICEOPS VERSION 1.0 LANGUAGES CXX)

# Set dylib paths relative to executable path under OSX
set(CMAKE_MACOSX_RPATH 0)
set(CMAKE_INSTALL_NAME_DIR @executable_path)
# Set attributes for Xcode project
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++17")
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_definitions(-D_THREAD_SAFE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
          "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
          FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# NOTE: do NOT add -DNDEBUG to Release.  The D original ran with asserts live
# in release builds, and sphere_database.cpp relies on assert(isFinite(...))
# guards; matching that behavior is part of the port contract.

# netCDF shared library (some systems like Homebrew install static and shared
# in the same place, so force the shared version)
find_library(NETCDF_LIBRARY NAMES netcdf
  PATHS /usr/local/lib/ /opt/local/lib/ /sw/lib/ /usr/lib64/ /usr/lib/ /usr/lib/x86_64-linux-gnu/)
find_path(NETCDF_INCLUDE_DIR NAMES netcdf.h
  PATHS /usr/local/include/ /opt/local/include/ /sw/include/ /usr/include/ /usr/include/x86_64-linux-gnu/)
if(NOT NETCDF_INCLUDE_DIR)
  message(FATAL_ERROR "netcdf.h not found; install the netCDF development headers.")
endif()

# CGAL is header-only since 5.0 and exports the CGAL::CGAL imported target.
# Don't include(${CGAL_USE_FILE}): UseCGAL.cmake was deprecated in CGAL 5 and
# removed in CGAL 6, where CGAL_USE_FILE is empty and the include() errors out.
find_package(CGAL REQUIRED COMPONENTS Core)
find_package(Eigen3 REQUIRED NO_MODULE)

# Third-party dependencies, fetched at configure time
include(FetchContent)

FetchContent_Declare(
  CLI11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG v2.4.2
)
FetchContent_Declare(
  nlohmann_json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.11.3
)
FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.5.4
)
# Python interpreter for the embedded matplotlib plotting facade
find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.13.6
)
FetchContent_MakeAvailable(CLI11 nlohmann_json Catch2 pybind11)

# Core library: noise database, trajectory, smoothing spline, plotting facade
add_library(niceops_core STATIC
  src/sphere_database.cpp
  src/trajectory.cpp
  src/plotting.cpp
  src/smoothing_spline.cpp
)
target_include_directories(niceops_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${NETCDF_INCLUDE_DIR}
)
target_link_libraries(niceops_core PUBLIC
  ${NETCDF_LIBRARY}
  Eigen3::Eigen
  nlohmann_json::nlohmann_json
  pybind11::embed
)

# CGAL 3D Delaunay interpolator over the operating-condition space
add_library(interpolator STATIC src/interpolator.cpp)
target_link_libraries(interpolator CGAL::CGAL)
target_link_libraries(niceops_core PUBLIC interpolator)

# Main executable
add_executable(niceops src/niceops.cpp)
target_link_libraries(niceops niceops_core CLI11::CLI11)

# Unit tests (Catch2)
add_executable(niceops_unit_tests
  tests/test_peak_range.cpp
  tests/test_sphere_database.cpp
  tests/test_trajectory.cpp
  tests/test_csv.cpp
  tests/test_legacy_repair.cpp
  tests/test_smoothing_spline.cpp
)
target_link_libraries(niceops_unit_tests niceops_core Catch2::Catch2WithMain)

# Golden output comparison tool (numeric, elementwise)
add_executable(golden_compare tests/golden_compare.cpp)

# Add testing with CTest
enable_testing()
# One default run feeds both the LAE and state golden comparisons
add_test(NAME regression_default COMMAND niceops -S state/ci_state.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# Validates every sphere in the database against its own stored levels
# (1.8M samples).  Was "broken_self_test" and excluded from CI while legacy
# hover directivity made it fail; it passes now and is a real test.
add_test(NAME self_test COMMAND niceops "--self_test" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME unit_test COMMAND niceops_unit_tests WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME regression_golden COMMAND golden_compare integrated/197421_LAE.csv golden/197421_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME state_golden COMMAND golden_compare state/ci_state.csv golden/state_197421.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
