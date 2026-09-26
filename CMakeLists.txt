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

# Link-time optimization.  The acoustic loop lives in niceops.cpp and calls
# Database::level_at_observer, which lives in sphere_database.cpp, ~100 million
# times a run; without LTO that is an out-of-line call across a translation-unit
# boundary every time, with the spills and the lost hoisting that implies.
# Turning it on is worth 19% of the loop and leaves the golden comparison
# passing unchanged, so it is on by default where the toolchain supports it.
option(NICEOPS_LTO "Build with link-time optimization" ON)
if(NICEOPS_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT NICEOPS_IPO_OK OUTPUT NICEOPS_IPO_WHY)
  if(NICEOPS_IPO_OK)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  else()
    message(STATUS "Link-time optimization unavailable: ${NICEOPS_IPO_WHY}")
  endif()
endif()

# netCDF shared library (some systems like Homebrew install static and shared
# in the same place, so force the shared version)
find_library(NETCDF_LIBRARY NAMES netcdf
  PATHS /usr/local/lib/ /opt/local/lib/ /sw/lib/ /usr/lib64/ /usr/lib/ /usr/lib/x86_64-linux-gnu/)
find_path(NETCDF_INCLUDE_DIR NAMES netcdf.h
  PATHS /usr/local/include/ /opt/local/include/ /sw/include/ /usr/include/ /usr/include/x86_64-linux-gnu/)
if(NOT NETCDF_INCLUDE_DIR)
  message(FATAL_ERROR "netcdf.h not found; install the netCDF development headers.")
endif()

# CGAL is gone: the Delaunay triangulation is now src/delaunay.cpp on the exact
# predicates in src/predicates.cpp, which also sheds the Boost, GMP and MPFR
# that CGAL pulled in behind it.

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
# Plotting is the only thing that needs Python.  Turning it off drops Python,
# pybind11 and matplotlib from the build entirely; the acoustics need none of
# them.
option(NICEOPS_PLOTTING "Build the matplotlib plotting facade" ON)
# The packages plotting needs are installed into a virtual environment under the
# build tree, so a run depends on what this build made rather than on whatever
# happens to be installed for the system interpreter.  Turn off to use the
# interpreter's own packages instead.
option(NICEOPS_BUNDLE_PYTHON "Create a virtual environment for matplotlib" ON)

if(NICEOPS_PLOTTING)
  find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11.git
    GIT_TAG v2.13.6
  )
  FetchContent_MakeAvailable(CLI11 nlohmann_json Catch2 pybind11)
else()
  FetchContent_MakeAvailable(CLI11 nlohmann_json Catch2)
endif()

# Build the virtual environment, once, and remember where its packages live.
# The embedded interpreter puts this directory at the front of sys.path.
set(NICEOPS_VENV_SITE_PACKAGES "")
if(NICEOPS_PLOTTING AND NICEOPS_BUNDLE_PYTHON)
  set(NICEOPS_VENV "${CMAKE_BINARY_DIR}/pyenv")
  if(WIN32)
    set(NICEOPS_VENV_PYTHON "${NICEOPS_VENV}/Scripts/python.exe")
  else()
    set(NICEOPS_VENV_PYTHON "${NICEOPS_VENV}/bin/python")
  endif()

  if(NOT EXISTS "${NICEOPS_VENV_PYTHON}")
    message(STATUS "Creating Python environment for plotting in ${NICEOPS_VENV}")
    execute_process(
      COMMAND "${Python3_EXECUTABLE}" -m venv "${NICEOPS_VENV}"
      RESULT_VARIABLE NICEOPS_VENV_RESULT
      OUTPUT_QUIET
    )
    if(NOT NICEOPS_VENV_RESULT EQUAL 0)
      message(WARNING "Could not create the plotting environment; falling back to the "
                      "system interpreter's packages.")
    endif()
  endif()

  if(EXISTS "${NICEOPS_VENV_PYTHON}")
    # Only reach the network when matplotlib is not already there, so repeat
    # configures stay offline and fast.
    execute_process(
      COMMAND "${NICEOPS_VENV_PYTHON}" -c "import matplotlib"
      RESULT_VARIABLE NICEOPS_HAVE_MATPLOTLIB
      OUTPUT_QUIET ERROR_QUIET
    )
    if(NOT NICEOPS_HAVE_MATPLOTLIB EQUAL 0)
      message(STATUS "Installing matplotlib into ${NICEOPS_VENV} (first configure only)")
      execute_process(
        COMMAND "${NICEOPS_VENV_PYTHON}" -m pip install --quiet --disable-pip-version-check
                matplotlib
        RESULT_VARIABLE NICEOPS_PIP_RESULT
      )
      if(NOT NICEOPS_PIP_RESULT EQUAL 0)
        message(WARNING "Could not install matplotlib. Plotting will report the failure at "
                        "run time; everything else is unaffected.")
      endif()
    endif()
    execute_process(
      COMMAND "${NICEOPS_VENV_PYTHON}" -c
              "import sysconfig; print(sysconfig.get_paths()['purelib'])"
      OUTPUT_VARIABLE NICEOPS_VENV_SITE_PACKAGES
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
    message(STATUS "Plotting packages: ${NICEOPS_VENV_SITE_PACKAGES}")
  endif()
endif()

if(NICEOPS_PLOTTING)
  set(NICEOPS_PLOTTING_SOURCE src/plotting.cpp)
else()
  set(NICEOPS_PLOTTING_SOURCE src/plotting_disabled.cpp)
endif()

# Core library: noise database, trajectory, smoothing spline, plotting facade
add_library(niceops_core STATIC
  src/sphere_database.cpp
  src/atmosphere.cpp
  src/weighting.cpp
  src/faddeeva.cpp
  src/ground_effect.cpp
  src/plate_receiver.cpp
  src/geodetic.cpp
  src/predicates.cpp
  src/delaunay.cpp
  # interpolator.cpp was its own library in the D era, to wall off CGAL's heavy
  # headers.  CGAL is gone, and the split had become a link cycle: interpolator
  # calls into Delaunay here, while niceops_core linked interpolator.  GNU ld
  # resolves archives in a single pass and failed on it; ld64 did not, so it
  # only ever broke Linux.
  src/interpolator.cpp
  src/trajectory.cpp
  ${NICEOPS_PLOTTING_SOURCE}
  src/smoothing_spline.cpp
)
target_include_directories(niceops_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${NETCDF_INCLUDE_DIR}
)
target_link_libraries(niceops_core PUBLIC
  ${NETCDF_LIBRARY}
  nlohmann_json::nlohmann_json
)
if(NICEOPS_PLOTTING)
  target_link_libraries(niceops_core PUBLIC pybind11::embed)
  if(NICEOPS_VENV_SITE_PACKAGES)
    set_source_files_properties(src/plotting.cpp PROPERTIES COMPILE_DEFINITIONS
      "NICEOPS_VENV_SITE_PACKAGES=\"${NICEOPS_VENV_SITE_PACKAGES}\"")
  endif()
endif()

# Main executable
add_executable(niceops src/niceops.cpp)
find_package(Threads REQUIRED)
target_link_libraries(niceops niceops_core CLI11::CLI11 Threads::Threads)

# Unit tests (Catch2)
add_executable(niceops_unit_tests
  tests/test_peak_range.cpp
  tests/test_sphere_database.cpp
  tests/test_geodetic.cpp
  tests/test_atmosphere.cpp
  tests/test_weighting.cpp
  tests/test_ground_effect.cpp
  tests/test_plate_receiver.cpp
  tests/test_spectral_propagation.cpp
  tests/test_predicates.cpp
  tests/test_delaunay.cpp
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
# The two golden comparisons read files this run writes.  CI happens to run the
# suite serially and in declaration order, but nothing said so: `ctest -j` was
# free to diff last run's output, or a file that did not exist yet.
set_tests_properties(regression_default PROPERTIES FIXTURES_SETUP demo_run)
# Validates every sphere in the database against its own stored levels
# (1.8M samples).  Was "broken_self_test" and excluded from CI while legacy
# hover directivity made it fail; it passes now and is a real test.
add_test(NAME self_test COMMAND niceops "--self_test" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME unit_test COMMAND niceops_unit_tests WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME regression_golden COMMAND golden_compare integrated/197421_LAE.csv golden/197421_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME state_golden COMMAND golden_compare state/ci_state.csv golden/state_197421.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(regression_golden state_golden PROPERTIES FIXTURES_REQUIRED demo_run)
# The spectral peak search and --db_cutoff skip what a bound puts below a
# threshold; with --ground_effect the bound includes the most the ground can
# add.  A skipped state must change nothing, so a run through the peak search
# has to equal a full march (-D forces one).  Coarse grid: this is about the
# wiring, which the bounds' own unit tests cannot see.
add_test(NAME ground_effect_windowed COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --ground_effect --turbulence 8e-4 --gridstep 492 -L warn -o integrated/ground_effect_windowed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ground_effect_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --ground_effect --turbulence 8e-4 --gridstep 492 -L warn -D integrated/ground_effect_dnl.csv -o integrated/ground_effect_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_effect_windowed ground_effect_marched PROPERTIES FIXTURES_SETUP ground_effect_runs)
add_test(NAME ground_effect_windowed_is_exact COMMAND golden_compare integrated/ground_effect_windowed.csv integrated/ground_effect_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_effect_windowed_is_exact PROPERTIES FIXTURES_REQUIRED ground_effect_runs)
# The same for the ground-plane receiver (--plate_table), whose bound is its
# own (PlateReceiver::gain_bounds), over the fitted 2017 Amedee ground.
add_test(NAME plate_windowed COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --gridstep 492 -L warn -o integrated/plate_windowed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME plate_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --gridstep 492 -L warn -D integrated/plate_dnl.csv -o integrated/plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(plate_windowed plate_marched PROPERTIES FIXTURES_SETUP plate_runs)
add_test(NAME plate_windowed_is_exact COMMAND golden_compare integrated/plate_windowed.csv integrated/plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(plate_windowed_is_exact PROPERTIES FIXTURES_REQUIRED plate_runs)
