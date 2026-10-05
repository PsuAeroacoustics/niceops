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

# assert() policy.  Shipped builds leave CMake's own -DNDEBUG in place, so the
# asserts in src/ compile out, as is usual for an optimized binary and about 4
# percent faster on the demo case.  Testing keeps them live: CI and the test
# recipe in CLAUDE.md pass -DNICEOPS_ASSERTS=ON, which strips NDEBUG from the
# optimized configurations and defines NICEOPS_ASSERTS_REQUIRED, which
# tests/test_build.cpp checks.  The asserts are side-effect free, so the two
# builds differ only where an assert would have fired.
option(NICEOPS_ASSERTS "Keep assert() live in optimized builds (for testing). Off for shipped binaries." OFF)
if (NICEOPS_ASSERTS)
  foreach(config RELEASE RELWITHDEBINFO MINSIZEREL)
    string(REPLACE "-DNDEBUG" "" CMAKE_CXX_FLAGS_${config} "${CMAKE_CXX_FLAGS_${config}}")
    string(STRIP "${CMAKE_CXX_FLAGS_${config}}" CMAKE_CXX_FLAGS_${config})
  endforeach()
  add_compile_definitions(NICEOPS_ASSERTS_REQUIRED)
endif()

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
  src/atmosphere_profile.cpp
  src/refraction.cpp
  src/stratified.cpp
  src/pade_propagator.cpp
  src/pe.cpp
  src/propagation_table.cpp
  src/airy.cpp
  src/exact_solutions.cpp
  src/weighting.cpp
  src/faddeeva.cpp
  src/ground_effect.cpp
  src/plate_receiver.cpp
  src/geodetic.cpp
  src/terrain.cpp
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
# The library itself starts threads (the observer loop, the PE's band pool and the
# propagation-table builds), so everything that links it needs the thread library;
# only niceops and niceops_benchmark said so, and the rest relied on a libc that
# folds pthread in.
find_package(Threads REQUIRED)
target_link_libraries(niceops_core PUBLIC
  ${NETCDF_LIBRARY}
  nlohmann_json::nlohmann_json
  Threads::Threads
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
target_link_libraries(niceops niceops_core CLI11::CLI11 Threads::Threads)

# Unit tests (Catch2)
add_executable(niceops_unit_tests
  tests/test_peak_range.cpp
  tests/test_parallel.cpp
  tests/test_build.cpp
  tests/test_iota.cpp
  tests/test_sphere_database.cpp
  tests/test_geodetic.cpp
  tests/test_terrain.cpp
  tests/test_atmosphere.cpp
  tests/test_atmosphere_profile.cpp
  tests/test_refraction.cpp
  tests/test_stratified.cpp
  tests/test_pade_propagator.cpp
  tests/test_propagation_table.cpp
  tests/test_pe.cpp
  tests/test_pe_shadow.cpp
  tests/test_pe_turbulence.cpp
  tests/test_benchmark.cpp
  tests/test_shadow_attenuation.cpp
  tests/test_shadow_transition.cpp
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
# Refraction models against the parabolic equation on one path (tools/pe_compare.cpp):
# a research tool, not part of a run.
add_executable(niceops_pe_compare tools/pe_compare.cpp)
target_link_libraries(niceops_pe_compare niceops_core CLI11::CLI11)
# The refracted ray's geometry for a batch of pairs, for evaluating ground models
# outside niceops on the same rays (tools/ray_geometry.cpp):
add_executable(niceops_ray_geometry tools/ray_geometry.cpp)
target_link_libraries(niceops_ray_geometry niceops_core CLI11::CLI11)
add_executable(niceops_pe_field tools/pe_field.cpp)
target_link_libraries(niceops_pe_field niceops_core CLI11::CLI11)
add_executable(niceops_propagation_table tools/propagation_table.cpp)
target_link_libraries(niceops_propagation_table niceops_core CLI11::CLI11)
# The Attenborough et al. (1995) benchmark cases by the PE and the exact
# solutions (tools/benchmark.cpp); tools/plot_benchmark.py draws them.
add_executable(niceops_benchmark tools/benchmark.cpp)
target_link_libraries(niceops_benchmark niceops_core CLI11::CLI11 Threads::Threads)

# Add testing with CTest
enable_testing()
add_test(NAME propagation_table_midpoint_validation COMMAND ${CMAKE_COMMAND}
  -DTOOL=$<TARGET_FILE:niceops_propagation_table>
  -DOUT=${CMAKE_CURRENT_BINARY_DIR}/table_validation
  -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_table_validation.cmake
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
# One default run feeds both the LAE and state golden comparisons
add_test(NAME regression_default COMMAND niceops -S state/ci_state.csv --no-spectral --no-ground_effect WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
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
# golden_compare is what every equivalence test below stands on, so it is tested
# itself.  A NaN or an infinity once compared equal to anything, which would have
# passed a footprint of nan.  The reject cases pass on the mismatch message
# (PASS_REGULAR_EXPRESSION replaces the exit status), so a usage error or a
# missing file cannot satisfy them.
add_test(NAME golden_compare_rejects_nan COMMAND golden_compare tests/data/compare_nan.csv tests/data/compare_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME golden_compare_rejects_inf COMMAND golden_compare tests/data/compare_inf.csv tests/data/compare_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME golden_compare_accepts_matching_inf COMMAND golden_compare tests/data/compare_same_inf.csv tests/data/compare_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(golden_compare_rejects_nan PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value 5: candidate nan")
set_tests_properties(golden_compare_rejects_inf PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value 8: candidate inf")
# -a mixes an ambient level in energy.  A negligible one (-200 dB) leaves the demo
# at the golden, which checks the plumbing of the marched, transformed-level path
# that -a forces; the arithmetic itself is a unit test (level_with_ambient).
add_test(NAME ambient_negligible_is_the_golden COMMAND niceops -S state/ci_state.csv --no-spectral --no-ground_effect -a -200 -o integrated/ambient_negligible.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ambient_negligible_matches COMMAND golden_compare integrated/ambient_negligible.csv golden/197421_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ambient_negligible_is_the_golden PROPERTIES FIXTURES_SETUP ambient_run)
set_tests_properties(ambient_negligible_matches PROPERTIES FIXTURES_REQUIRED ambient_run)
# A failure used to reach only log/niceops.log, so a skipped trajectory left the
# terminal empty with exit status 1 as the only sign.  Warnings and worse are now
# echoed to stderr; the case passes on that text.
add_test(NAME missing_track_is_reported COMMAND niceops -t tracks/does_not_exist.csv -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(missing_track_is_reported PROPERTIES PASS_REGULAR_EXPRESSION "Skipping tracks/does_not_exist.csv: Unable to open CSV file")
# --export_sphere writes with NC_CLOBBER: giving it the database it reads replaced that
# database with a one-sphere file and exited 0.
add_test(NAME export_guard_copy COMMAND ${CMAKE_COMMAND} -E copy databases/S-76D_fixed_load_factor.nod integrated/export_guard.nod WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME export_over_the_database_is_refused COMMAND niceops -d integrated/export_guard.nod --export_sphere integrated/export_guard.nod --sphere_condition 0.1 0 0.006 -S state/ci_state.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(export_guard_copy PROPERTIES FIXTURES_SETUP export_guard_file)
set_tests_properties(export_over_the_database_is_refused PROPERTIES FIXTURES_REQUIRED export_guard_file PASS_REGULAR_EXPRESSION "would overwrite the database it reads")
# A negative --receiver_height was taken as typed (the ground model assumes the microphone is
# at or above the plane, and its table disagrees with the exact model below it).
add_test(NAME negative_receiver_height_is_refused COMMAND niceops --ground_effect --receiver_height -3.937 -S state/ci_state.csv -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(negative_receiver_height_is_refused PROPERTIES PASS_REGULAR_EXPRESSION "--receiver_height: Value -3.937 not in range")
# An exception in an observer-loop worker must skip the trajectory, not abort the
# process.  Making the first history file's name a directory makes opening it
# throw inside a worker; the case passes only on the "Skipping" line that
# run()'s handler writes, which a std::terminate never reaches.
add_test(NAME worker_exception_clean_clean COMMAND ${CMAKE_COMMAND} -E remove_directory integrated/blocked_hist WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME worker_exception_clean_block COMMAND ${CMAKE_COMMAND} -E make_directory integrated/blocked_hist/h_i0.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME worker_exception_is_skipped COMMAND niceops -S state/ci_state.csv --gridstep 1640 -j 4 -i integrated/blocked_hist/h -o integrated/blocked_hist_levels.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(worker_exception_clean_clean PROPERTIES FIXTURES_SETUP blocked_hist_clean)
set_tests_properties(worker_exception_clean_block PROPERTIES FIXTURES_REQUIRED blocked_hist_clean FIXTURES_SETUP blocked_hist_ready)
set_tests_properties(worker_exception_is_skipped PROPERTIES FIXTURES_REQUIRED blocked_hist_ready PASS_REGULAR_EXPRESSION "Skipping tracks/S-76D/197421AC.csv: Unable to open file for writing")
# With several trajectories every observer index restarts at 0, so -i and
# --band_history files named only by index were overwritten by each trajectory in
# turn.  They are now named for their track as well; and two trajectories asked
# to write one file are refused up front rather than silently overwriting.
add_test(NAME multi_trajectory_history_clean COMMAND ${CMAKE_COMMAND} -E remove_directory integrated/multi_hist WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME multi_trajectory_history_run COMMAND niceops -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197405AC.csv --gridstep 1640 -L warn -i integrated/multi_hist/h -o integrated/multi_a.csv -o integrated/multi_b.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME multi_trajectory_history_names COMMAND ${CMAKE_COMMAND} "-DFILES=integrated/multi_hist/h_197421AC_i0.csv;integrated/multi_hist/h_197405AC_i0.csv;integrated/multi_hist/h_197421AC_i48.csv;integrated/multi_hist/h_197405AC_i48.csv" "-DABSENT=integrated/multi_hist/h_i0.csv" -P tests/check_files.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(multi_trajectory_history_clean PROPERTIES FIXTURES_SETUP multi_history_clean)
set_tests_properties(multi_trajectory_history_run PROPERTIES FIXTURES_REQUIRED multi_history_clean FIXTURES_SETUP multi_history_run)
set_tests_properties(multi_trajectory_history_names PROPERTIES FIXTURES_REQUIRED multi_history_run)
add_test(NAME duplicate_output_names_rejected COMMAND niceops -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197405AC.csv --gridstep 1640 -o integrated/same.csv -o integrated/same.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(duplicate_output_names_rejected PROPERTIES PASS_REGULAR_EXPRESSION "Two trajectories would write the same output file")
# -D on a geodetic track with no utcsec column and no --start_utc used to date the
# flight at UTC midnight and add the night penalty (+10 dB) without a word.  It is
# refused now, and a start time makes the same run work.
add_test(NAME dnl_needs_a_clock COMMAND niceops -t tests/data/demo_geodetic_noclock.csv --origin 39.82914323 -74.40996141 --frame_bearing 210 --ground_elevation 54.2651 --gridstep 1640 -D integrated/noclock_dnl.csv -o integrated/noclock_levels.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(dnl_needs_a_clock PROPERTIES PASS_REGULAR_EXPRESSION "-D needs the time of day")
add_test(NAME dnl_with_start_utc COMMAND niceops -t tests/data/demo_geodetic_noclock.csv --origin 39.82914323 -74.40996141 --frame_bearing 210 --ground_elevation 54.2651 --gridstep 1640 --start_utc 43200 -D integrated/dated_dnl.csv -o integrated/dated_levels.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# An observer exactly at a source position: spherical spreading diverges at r = 0.
# The level was +inf in the CSV, and with the asserts live the run aborted.  The
# case passes only if the run finishes.
add_test(NAME observer_on_the_source_survives COMMAND niceops -b tests/data/observer_on_source.csv -o integrated/observer_on_source.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# The vehicle file against the database it is run with.  Only the tip speed was
# ever compared, in an info line; weight and radius set CT in the same way, and a
# vehicle file 10 percent heavy moved the footprint 0.75 dB median with nothing
# said.  All three are warnings now.  The B407 file and database are a matched
# pair (the weights differ by panam's g = 9.82, 0.14 percent) and must stay quiet.
set(B407_RUN -d databases/B407_fixed_load_factor.nod -t tracks/S-76D/197421AC.csv --gridstep 1640 --humidity 20 -o integrated/vehicle_check.csv)
add_test(NAME vehicle_tip_speed_mismatch_warns COMMAND niceops ${B407_RUN} -v databases/S-76D.json --accept_vehicle_mismatch WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(vehicle_tip_speed_mismatch_warns PROPERTIES PASS_REGULAR_EXPRESSION "Vehicle main rotor tip speed disagrees")
add_test(NAME vehicle_weight_mismatch_warns COMMAND niceops ${B407_RUN} -v tests/data/B407_heavy.json --accept_vehicle_mismatch WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME vehicle_weight_mismatch_refused COMMAND niceops ${B407_RUN} -v tests/data/B407_heavy.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(vehicle_weight_mismatch_refused PROPERTIES PASS_REGULAR_EXPRESSION "percent off")
add_test(NAME lamax_file_is_written COMMAND sh -c "rm -f integrated/lamax.csv && $<TARGET_FILE:niceops> -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --gridstep 1640 --lamax_file integrated/lamax.csv -o integrated/lamax_run.csv && test -s integrated/lamax.csv" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME static_hover_runs COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tests/data/static_hover.csv --gridstep 1640 -L warn -o integrated/static_hover.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME nax_file_is_written COMMAND sh -c "rm -f integrated/nax.csv && $<TARGET_FILE:niceops> -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --gridstep 1640 --nax_file integrated/nax.csv --nax_threshold 90 -o integrated/nax_run.csv && test -s integrated/nax.csv" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME nax_needs_a_threshold COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --gridstep 1640 --nax_file integrated/nax_unset.csv -o integrated/nax_unset_run.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nax_needs_a_threshold PROPERTIES WILL_FAIL TRUE)
add_test(NAME full_history_runs COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --gridstep 1640 --full_history -o integrated/full_history.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME load_factor_is_held_at_the_ceiling COMMAND niceops ${B407_RUN} -v databases/B407.json --max_load_factor 1 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(load_factor_is_held_at_the_ceiling PROPERTIES PASS_REGULAR_EXPRESSION "exceed 1 g on a fixed_load_factor database")
set_tests_properties(vehicle_weight_mismatch_warns PROPERTIES PASS_REGULAR_EXPRESSION "Vehicle weight disagrees")
add_test(NAME matched_vehicle_is_silent COMMAND niceops ${B407_RUN} -v databases/B407.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(matched_vehicle_is_silent PROPERTIES FAIL_REGULAR_EXPRESSION "disagrees;Relative humidity is not set")
# Humidity is the one atmosphere value that is usually defaulted (the vehicle
# files give a temperature and a density, not a humidity) and absorption is very
# sensitive to it, so a defaulted one is said aloud; the other values say where
# they came from in the log.
add_test(NAME humidity_default_is_warned COMMAND niceops -d databases/B407_fixed_load_factor.nod -t tracks/S-76D/197421AC.csv --gridstep 1640 -v databases/B407.json -o integrated/humidity_default.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(humidity_default_is_warned PROPERTIES PASS_REGULAR_EXPRESSION "Relative humidity is not set")
# Unit slips in the atmosphere options are refused before anything loads.
add_test(NAME temperature_in_celsius_rejected COMMAND niceops --temperature 15 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(temperature_in_celsius_rejected PROPERTIES PASS_REGULAR_EXPRESSION "--temperature: a temperature of 15.* is it in Celsius")
add_test(NAME pressure_in_pascals_rejected COMMAND niceops --pressure 101325 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(pressure_in_pascals_rejected PROPERTIES PASS_REGULAR_EXPRESSION "--pressure: a pressure of 101325.* is it in Pa or hPa")
add_test(NAME humidity_as_fraction_rejected COMMAND niceops --humidity 0.5 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(humidity_as_fraction_rejected PROPERTIES PASS_REGULAR_EXPRESSION "--humidity: a relative humidity of 0.5.* not as a fraction")
# An observer file with a row that does not parse, or with rows of different
# widths, used to run on the rows that did: a typo removed an observer, a short
# row became an observer at the wrong coordinates.
add_test(NAME observer_typo_rejected COMMAND niceops -b tests/data/observers_typo.csv -o integrated/observers_typo.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(observer_typo_rejected PROPERTIES PASS_REGULAR_EXPRESSION "line 3 is not a row of numbers")
add_test(NAME observer_mixed_widths_rejected COMMAND niceops -b tests/data/observers_mixed_widths.csv -o integrated/observers_mixed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(observer_mixed_widths_rejected PROPERTIES PASS_REGULAR_EXPRESSION "has 2 fields where the first has 3")
# The PE tools took numeric options unchecked: --source_height 0 without --top
# looped until it had eaten 21 GB, --max_range 0 and --range_step 0 threw inside a
# worker thread (an abort), --samples 0 divided by zero and wrote a table of
# +Infinity with exit status 0, and --max_angle above the documented 60 degrees
# reported paths the march cannot carry instead of NaN.  Each is refused.
set(PE_FIELD_BASE --atmosphere atmospheres/example_windy.csv --toward 312)
add_test(NAME pe_field_zero_source_height_needs_top COMMAND niceops_pe_field ${PE_FIELD_BASE} --source_height 0 --bands 500 --max_range 2000 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(pe_field_zero_source_height_needs_top PROPERTIES PASS_REGULAR_EXPRESSION "gives no default for --top")
add_test(NAME pe_field_rejects_zero_range COMMAND niceops_pe_field ${PE_FIELD_BASE} --source_height 100 --max_range 0 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(pe_field_rejects_zero_range PROPERTIES PASS_REGULAR_EXPRESSION "--max_range: Value 0 not in range")
add_test(NAME pe_field_rejects_zero_samples COMMAND niceops_pe_field ${PE_FIELD_BASE} --source_height 100 --max_range 1000 --samples 0 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(pe_field_rejects_zero_samples PROPERTIES PASS_REGULAR_EXPRESSION "--samples: Value 0 not in range")
add_test(NAME pe_field_rejects_steep_max_angle COMMAND niceops_pe_field ${PE_FIELD_BASE} --source_height 100 --max_range 1000 --max_angle 75 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(pe_field_rejects_steep_max_angle PROPERTIES PASS_REGULAR_EXPRESSION "--max_angle: Value 75 not in range")
# The spectral peak search and --db_cutoff skip what a bound puts below a
# threshold; with --ground_effect the bound includes the most the ground can
# add.  A skipped state must change nothing, so a run through the peak search
# has to equal a full march (-D forces one).  Coarse grid: this is about the
# wiring, which the bounds' own unit tests cannot see.
add_test(NAME ground_effect_windowed COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --ground_effect --turbulence 8e-4 --gridstep 492 -L warn -o integrated/ground_effect_windowed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ground_effect_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --ground_effect --turbulence 8e-4 --gridstep 492 -L warn -D integrated/ground_effect_dnl.csv -o integrated/ground_effect_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_effect_windowed ground_effect_marched PROPERTIES FIXTURES_SETUP ground_effect_runs)
add_test(NAME ground_effect_windowed_is_exact COMMAND golden_compare integrated/ground_effect_windowed.csv integrated/ground_effect_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_effect_windowed_is_exact PROPERTIES FIXTURES_REQUIRED ground_effect_runs)
# The same for the ground-plane receiver (--plate_table), whose bound is its
# own (PlateReceiver::gain_bounds), over the fitted 2017 Amedee ground.
add_test(NAME plate_windowed COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --gridstep 492 -L warn -o integrated/plate_windowed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME plate_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --gridstep 492 -L warn -D integrated/plate_dnl.csv -o integrated/plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(plate_windowed plate_marched PROPERTIES FIXTURES_SETUP plate_runs)
add_test(NAME plate_windowed_is_exact COMMAND golden_compare integrated/plate_windowed.csv integrated/plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(plate_windowed_is_exact PROPERTIES FIXTURES_REQUIRED plate_runs)
# --propagation_table, end to end: a table built by niceops_propagation_table
# from the straight-ray ground model must give the footprint --ground_effect
# gives, to within the table's interpolation (10 ft in height and range: 0.19 dB
# at worst on this grid, 0.02 at 5 ft).  It checks the generator, the file, the lookup and the
# footprint's use of it together.  Source heights span the demo track's.
add_test(NAME propagation_table_build COMMAND niceops_propagation_table --method ground_model
    --database databases/B407_fixed_load_factor.nod --sound_speed 335.636 --flow_resistance 225
    --source_heights 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350 360 370 380 390 400 410 420 430 440 450 460 470 480 490 500 510 520 530 540 550 560 570 580 590 600 610 620 630 640 650 660 670 680 690 700 710 720 730 740 750 760 770 780 790 800 810 820 830 840 850 860 870 880 890 900 910 920 930 940 950 960 970 980 990 1000 1010 1020 1030 1040 1050 1060 1070 1080 1090 1100 1110 1120 1130
    --receiver_heights 3.93701 --max_range 15000 --step 10 -o tables/ground_model_test.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_run COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --propagation_table tables/ground_model_test.nc --gridstep 492 -L warn -o integrated/propagation_table.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_reference COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --ground_effect --exact_ground --gridstep 492 -L warn -o integrated/propagation_table_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_build PROPERTIES FIXTURES_SETUP propagation_table_file)
set_tests_properties(propagation_table_run PROPERTIES FIXTURES_REQUIRED propagation_table_file FIXTURES_SETUP propagation_table_runs)
set_tests_properties(propagation_table_reference PROPERTIES FIXTURES_SETUP propagation_table_runs)
# The peak search with a table: its bound takes the table's largest gain per
# band, so a windowed run must equal a marched one (-D forces the march) to the bit.
add_test(NAME propagation_table_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --propagation_table tables/ground_model_test.nc --gridstep 492 -L warn -D integrated/propagation_table_dnl.csv -o integrated/propagation_table_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_marched PROPERTIES FIXTURES_REQUIRED propagation_table_file FIXTURES_SETUP propagation_table_runs)
add_test(NAME propagation_table_windowed_is_exact COMMAND golden_compare integrated/propagation_table.csv integrated/propagation_table_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_windowed_is_exact PROPERTIES FIXTURES_REQUIRED propagation_table_runs)
# A table with a +25 dB bump beyond 1500 m and a bearing tilt, on eight bearings, run with the frame
# turned 210 degrees: the peak search's bound has to carry the table's largest gain, so the windowed run
# equals the marched one; and a table that says it carries absorption makes the run add none, so the
# footprint is the louder of the two by the absorption (golden_compare reports the mismatch).
add_test(NAME propagation_table_bump_build COMMAND niceops_unit_tests "[.bump_table]" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_bump_build PROPERTIES FIXTURES_REQUIRED propagation_table_file FIXTURES_SETUP propagation_table_bump_file)
add_test(NAME propagation_table_bump COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --propagation_table tables/bump_test.nc --frame_bearing 210 --gridstep 492 -L warn -o integrated/propagation_table_bump.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_bump_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --propagation_table tables/bump_test.nc --frame_bearing 210 --gridstep 492 -L warn -D integrated/propagation_table_bump_dnl.csv -o integrated/propagation_table_bump_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_bump_absorbed COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --propagation_table tables/bump_absorbed_test.nc --frame_bearing 210 --gridstep 492 -L warn -o integrated/propagation_table_bump_absorbed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_bump propagation_table_bump_marched propagation_table_bump_absorbed PROPERTIES FIXTURES_REQUIRED propagation_table_bump_file FIXTURES_SETUP propagation_table_bump_runs)
add_test(NAME propagation_table_bump_windowed_is_exact COMMAND golden_compare integrated/propagation_table_bump.csv integrated/propagation_table_bump_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_absorption_is_not_added_twice COMMAND golden_compare integrated/propagation_table_bump.csv integrated/propagation_table_bump_absorbed.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_bump_windowed_is_exact PROPERTIES FIXTURES_REQUIRED propagation_table_bump_runs)
set_tests_properties(propagation_table_absorption_is_not_added_twice PROPERTIES FIXTURES_REQUIRED propagation_table_bump_runs PASS_REGULAR_EXPRESSION "mismatch at value")
add_test(NAME propagation_table_is_the_ground_model COMMAND golden_compare integrated/propagation_table.csv integrated/propagation_table_reference.csv 0.25 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_is_the_ground_model PROPERTIES FIXTURES_REQUIRED propagation_table_runs)
# --atmosphere through a single-level profile holding the vehicle file's own
# atmosphere has to reproduce the homogeneous run: the whole layered-absorption
# path (the table, the per-state and per-observer integrals, the per-path mean)
# against the one run-wide coefficient array.  Not bit-identical -- a mean of a
# uniform integral is the coefficient to a few ULP -- so golden_compare's 1e-9.
add_test(NAME atmosphere_uniform COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/B407_uniform.csv --gridstep 492 -L warn -o integrated/atmosphere_uniform.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME atmosphere_homogeneous COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --gridstep 492 -L warn -o integrated/atmosphere_homogeneous.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(atmosphere_uniform atmosphere_homogeneous PROPERTIES FIXTURES_SETUP atmosphere_runs)
add_test(NAME atmosphere_uniform_is_homogeneous COMMAND golden_compare integrated/atmosphere_uniform.csv integrated/atmosphere_homogeneous.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(atmosphere_uniform_is_homogeneous PROPERTIES FIXTURES_REQUIRED atmosphere_runs)
# --speed_reference air through the same atmosphere with its wind given as zero:
# every flight condition is the ground one, so it has to reproduce the uniform
# run.  And an air reference with no wind to subtract is refused, not taken as
# calm.
add_test(NAME speed_reference_air_calm COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere tests/data/B407_uniform_calm.csv --speed_reference air --gridstep 492 -L warn -o integrated/speed_reference_air_calm.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(speed_reference_air_calm PROPERTIES FIXTURES_SETUP atmosphere_runs)
add_test(NAME speed_reference_air_calm_is_ground COMMAND golden_compare integrated/speed_reference_air_calm.csv integrated/atmosphere_uniform.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(speed_reference_air_calm_is_ground PROPERTIES FIXTURES_REQUIRED atmosphere_runs)
add_test(NAME speed_reference_air_needs_wind COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/B407_uniform.csv --speed_reference air --gridstep 492 -L warn -o integrated/speed_reference_air_needs_wind.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(speed_reference_air_needs_wind PROPERTIES WILL_FAIL TRUE)
# --nose.  The hidden Catch case [.heading_tracks] writes two copies of the demo
# track differing only in the heading column -- the ground track's own compass
# bearing, and that turned 15 degrees clockwise -- and a one-sphere database whose
# root says azimuth_reference = heading.  Heading mode along the track has to
# reproduce the default nose; crabbed, it has to move the footprint; and the
# database's attribute has to select heading mode by itself.
add_test(NAME nose_heading_build COMMAND niceops_unit_tests "[.heading_tracks]" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nose_heading_build PROPERTIES FIXTURES_SETUP heading_tracks)
add_test(NAME nose_default_on_track COMMAND niceops -t integrated/heading_on_track.csv --frame_bearing 210 --gridstep 328 -L warn -o integrated/nose_default_on_track.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME nose_heading_on_track COMMAND niceops -t integrated/heading_on_track.csv --frame_bearing 210 --gridstep 328 --nose heading -L warn -o integrated/nose_heading_on_track.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME nose_heading_crabbed COMMAND niceops -t integrated/heading_crabbed.csv --frame_bearing 210 --gridstep 328 --nose heading -L warn -o integrated/nose_heading_crabbed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME nose_database_heading COMMAND niceops -d integrated/heading_sphere.nod -t integrated/heading_crabbed.csv --frame_bearing 210 --gridstep 328 -L warn -o integrated/nose_database_heading.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME nose_database_heading_explicit COMMAND niceops -d integrated/heading_sphere.nod -t integrated/heading_crabbed.csv --frame_bearing 210 --gridstep 328 --nose heading -L warn -o integrated/nose_database_heading_explicit.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME nose_database_velocity COMMAND niceops -d integrated/heading_sphere.nod -t integrated/heading_crabbed.csv --frame_bearing 210 --gridstep 328 --nose velocity -L warn -o integrated/nose_database_velocity.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nose_database_velocity PROPERTIES PASS_REGULAR_EXPRESSION "--nose velocity overrides the database, whose sphere azimuth is measured from the heading")
set_tests_properties(nose_default_on_track nose_heading_on_track nose_heading_crabbed nose_database_heading nose_database_heading_explicit nose_database_velocity PROPERTIES FIXTURES_REQUIRED heading_tracks FIXTURES_SETUP nose_runs)
add_test(NAME nose_heading_on_track_is_velocity COMMAND golden_compare integrated/nose_heading_on_track.csv integrated/nose_default_on_track.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# A crab of 15 degrees moves the footprint (golden_compare reports a mismatch, which is the pass).
add_test(NAME nose_heading_crabbed_differs COMMAND golden_compare integrated/nose_heading_crabbed.csv integrated/nose_default_on_track.csv 0.1 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nose_heading_crabbed_differs PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value")
add_test(NAME nose_database_heading_is_default COMMAND golden_compare integrated/nose_database_heading.csv integrated/nose_database_heading_explicit.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME nose_database_velocity_differs COMMAND golden_compare integrated/nose_database_velocity.csv integrated/nose_database_heading.csv 0.1 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nose_database_velocity_differs PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value")
set_tests_properties(nose_heading_on_track_is_velocity nose_heading_crabbed_differs nose_database_heading_is_default nose_database_velocity_differs PROPERTIES FIXTURES_REQUIRED nose_runs)
# The hover nose: with --frame_bearing given, two tracks whose hovers share a
# heading and whose approaches differ only in their crab (which the default nose
# ignores, flying along the velocity) give the same footprint.  The track's own
# fit, which used to win over the frame, absorbed the crabbed approach's 15
# degrees and turned that hover's nose by them.
add_test(NAME nose_hover_crabbed COMMAND niceops -t integrated/heading_hover_crabbed.csv --frame_bearing 210 --gridstep 328 -L warn -o integrated/nose_hover_crabbed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME nose_hover_straight COMMAND niceops -t integrated/heading_hover_straight.csv --frame_bearing 210 --gridstep 328 -L warn -o integrated/nose_hover_straight.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nose_hover_crabbed nose_hover_straight PROPERTIES FIXTURES_REQUIRED heading_tracks FIXTURES_SETUP nose_hover_runs)
add_test(NAME nose_hover_takes_the_frame COMMAND golden_compare integrated/nose_hover_crabbed.csv integrated/nose_hover_straight.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nose_hover_takes_the_frame PROPERTIES FIXTURES_REQUIRED nose_hover_runs)
# Heading mode is refused with no frame to turn the heading by (the track's own
# fit would absorb the crab), with a database implying it as with --nose, and on
# a track with no heading at all.
add_test(NAME nose_heading_needs_a_frame COMMAND niceops -t integrated/heading_crabbed.csv --nose heading --gridstep 1640 -L warn -o integrated/nose_heading_no_frame.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nose_heading_needs_a_frame PROPERTIES PASS_REGULAR_EXPRESSION "--nose heading turns each state's compass heading into the local frame, and needs that frame's bearing" FIXTURES_REQUIRED heading_tracks)
add_test(NAME nose_database_heading_needs_a_frame COMMAND niceops -d integrated/heading_sphere.nod -t integrated/heading_crabbed.csv --gridstep 1640 -L warn -o integrated/nose_database_no_frame.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nose_database_heading_needs_a_frame PROPERTIES PASS_REGULAR_EXPRESSION "azimuth_reference \\(heading\\) turns each state's compass heading.*or --nose velocity" FIXTURES_REQUIRED heading_tracks)
add_test(NAME nose_heading_needs_a_heading COMMAND niceops -t tests/data/demo_geodetic_noclock.csv --origin 39.82914323 -74.40996141 --frame_bearing 210 --ground_elevation 54.2651 --gridstep 1640 --nose heading -L warn -o integrated/nose_heading_no_heading.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nose_heading_needs_a_heading PROPERTIES PASS_REGULAR_EXPRESSION "--nose heading needs the nose's compass heading, and tests/data/demo_geodetic_noclock.csv carries none")
# A database implying heading mode on a track with no heading falls back to the
# velocity, with a warning, and runs.
add_test(NAME nose_database_heading_without_heading COMMAND niceops -d integrated/heading_sphere.nod -t tests/data/demo_geodetic_noclock.csv --origin 39.82914323 -74.40996141 --frame_bearing 210 --ground_elevation 54.2651 --gridstep 1640 -L warn -o integrated/nose_database_no_heading.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nose_database_heading_without_heading PROPERTIES PASS_REGULAR_EXPRESSION "the nose follows the velocity instead" FIXTURES_REQUIRED heading_tracks)
# On the spectral path the windowed peak search bounds each state through its
# sphere frame, so a heading nose has to leave it exact.
add_test(NAME nose_heading_spectral COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t integrated/heading_crabbed.csv --frame_bearing 210 --nose heading --gridstep 984 -L warn -o integrated/nose_heading_spectral.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME nose_heading_spectral_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t integrated/heading_crabbed.csv --frame_bearing 210 --nose heading --gridstep 984 -L warn -D integrated/nose_heading_spectral_dnl.csv -o integrated/nose_heading_spectral_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nose_heading_spectral nose_heading_spectral_marched PROPERTIES FIXTURES_REQUIRED heading_tracks FIXTURES_SETUP nose_spectral_runs)
add_test(NAME nose_heading_windowed_is_exact COMMAND golden_compare integrated/nose_heading_spectral.csv integrated/nose_heading_spectral_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nose_heading_windowed_is_exact PROPERTIES FIXTURES_REQUIRED nose_spectral_runs)
# --nose ground_velocity keeps the nose on the ground track under --speed_reference
# air, which otherwise turns it onto the air velocity: in a wind the two differ.
add_test(NAME nose_air COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --frame_bearing 210 --speed_reference air --gridstep 984 -L warn -o integrated/nose_air.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME nose_air_ground_velocity COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --frame_bearing 210 --speed_reference air --nose ground_velocity --gridstep 984 -L warn -o integrated/nose_air_ground_velocity.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nose_air nose_air_ground_velocity PROPERTIES FIXTURES_SETUP nose_air_runs)
add_test(NAME nose_ground_velocity_is_not_air COMMAND golden_compare integrated/nose_air_ground_velocity.csv integrated/nose_air.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(nose_ground_velocity_is_not_air PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value" FIXTURES_REQUIRED nose_air_runs)
# The same over arrays/coyle_approach.csv, whose microphones stand at different
# heights.  A grid's observers share one, so their mean coefficients are formed
# once per state; these take the per-pair path the grid never reaches.
add_test(NAME atmosphere_uniform_scattered COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv -b arrays/coyle_approach.csv --atmosphere atmospheres/B407_uniform.csv -L warn -o integrated/atmosphere_uniform_scattered.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME atmosphere_homogeneous_scattered COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -b arrays/coyle_approach.csv -L warn -o integrated/atmosphere_homogeneous_scattered.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(atmosphere_uniform_scattered atmosphere_homogeneous_scattered PROPERTIES FIXTURES_SETUP atmosphere_scattered_runs)
add_test(NAME atmosphere_uniform_scattered_is_homogeneous COMMAND golden_compare integrated/atmosphere_uniform_scattered.csv integrated/atmosphere_homogeneous_scattered.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(atmosphere_uniform_scattered_is_homogeneous PROPERTIES FIXTURES_REQUIRED atmosphere_scattered_runs)
# --refraction linear through the uniform profile: no gradient, so every path is
# straight, and it has to reproduce the homogeneous run above -- the whole
# refracted-path plumbing (per-path fits, travel time into reception time, the
# march it forces) against the straight path.
add_test(NAME refraction_uniform COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/B407_uniform.csv --refraction linear --gridstep 492 -L warn -o integrated/refraction_uniform.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(refraction_uniform PROPERTIES FIXTURES_SETUP atmosphere_runs)
add_test(NAME refraction_uniform_is_homogeneous COMMAND golden_compare integrated/refraction_uniform.csv integrated/atmosphere_homogeneous.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(refraction_uniform_is_homogeneous PROPERTIES FIXTURES_REQUIRED atmosphere_runs)
# The configuration most test data is compared in: flush plate microphones,
# scattered over sloping ground (arrays/coyle_approach.csv stands above and
# below the plane).  Refraction in uniform air must end each ray where the
# straight model does.
add_test(NAME refraction_uniform_plate_scattered COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv -b arrays/coyle_approach.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --atmosphere atmospheres/B407_uniform.csv --refraction linear -L warn -o integrated/refraction_uniform_plate_scattered.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME straight_plate_scattered COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -b arrays/coyle_approach.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc -L warn -o integrated/straight_plate_scattered.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(refraction_uniform_plate_scattered straight_plate_scattered PROPERTIES FIXTURES_SETUP refraction_plate_runs)
add_test(NAME refraction_uniform_plate_is_straight COMMAND golden_compare integrated/refraction_uniform_plate_scattered.csv integrated/straight_plate_scattered.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(refraction_uniform_plate_is_straight PROPERTIES FIXTURES_REQUIRED refraction_plate_runs)
# One run over an array of flush plate microphones and microphones on poles at
# two heights (id,x,y,z,height) has to equal a run over each kind alone: the
# plates with --plate_table, each pole height as the run's --receiver_height.
# --exact_ground, because a table's extent follows its run's observers.  Then
# the mixed run through the peak search against a march (-D), since its bound
# is the larger of the plate's and the point receiver's.
add_test(NAME mixed_receivers COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --exact_ground -L warn --plate_table plates/gr1425_vp200.nc -b tests/data/mixed_receivers.csv -o integrated/mixed_receivers.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME mixed_receivers_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --exact_ground -L warn --plate_table plates/gr1425_vp200.nc -b tests/data/mixed_receivers.csv -D integrated/mixed_receivers_dnl.csv -o integrated/mixed_receivers_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME mixed_flush_only COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --exact_ground -L warn --plate_table plates/gr1425_vp200.nc -b tests/data/mixed_flush.csv -o integrated/mixed_flush.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME mixed_poles_4ft_only COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --exact_ground -L warn --receiver_height 3.937 -b tests/data/mixed_poles_4ft.csv -o integrated/mixed_poles_4ft.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME mixed_poles_10ft_only COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --exact_ground -L warn --receiver_height 10 -b tests/data/mixed_poles_10ft.csv -o integrated/mixed_poles_10ft.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(mixed_receivers mixed_receivers_marched mixed_flush_only mixed_poles_4ft_only mixed_poles_10ft_only PROPERTIES FIXTURES_SETUP mixed_runs)
add_test(NAME mixed_expected COMMAND ${CMAKE_COMMAND} -DOUT=integrated/mixed_expected.csv "-DINPUTS=integrated/mixed_flush.csv;integrated/mixed_poles_4ft.csv;integrated/mixed_poles_10ft.csv" -P tests/concatenate.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(mixed_expected PROPERTIES FIXTURES_REQUIRED mixed_runs FIXTURES_SETUP mixed_expected_file)
add_test(NAME mixed_receivers_are_their_subsets COMMAND golden_compare integrated/mixed_receivers.csv integrated/mixed_expected.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(mixed_receivers_are_their_subsets PROPERTIES FIXTURES_REQUIRED "mixed_runs;mixed_expected_file")
add_test(NAME mixed_receivers_windowed_is_exact COMMAND golden_compare integrated/mixed_receivers.csv integrated/mixed_receivers_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(mixed_receivers_windowed_is_exact PROPERTIES FIXTURES_REQUIRED mixed_runs)
# --band_history over flush plate microphones.  Writing the bands forces a
# march, and must leave the footprint equal to the peak search's; a refracted
# path through uniform air must record the straight one's geometry and bands.
# --band_limits spans every band, so no field is left empty for golden_compare
# to skip the row over.  (That the bands re-weight to the level is a unit test,
# "band_levels are the bands the weighted level sums".)
add_test(NAME band_history_straight COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -b tests/data/mixed_flush.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --band_limits 5 20000 -L warn --band_history integrated/band_history/straight -o integrated/band_history_straight.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME band_history_refracted COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv -b tests/data/mixed_flush.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --band_limits 5 20000 -L warn --atmosphere atmospheres/B407_uniform.csv --refraction linear --band_history integrated/band_history/refracted -o integrated/band_history_refracted.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME band_history_windowed COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -b tests/data/mixed_flush.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --band_limits 5 20000 -L warn -o integrated/band_history_windowed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(band_history_straight band_history_refracted band_history_windowed PROPERTIES FIXTURES_SETUP band_history_runs)
add_test(NAME band_history_leaves_footprint COMMAND golden_compare integrated/band_history_straight.csv integrated/band_history_windowed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME band_history_refracted_is_straight COMMAND golden_compare integrated/band_history/refracted_i3.csv integrated/band_history/straight_i3.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(band_history_leaves_footprint band_history_refracted_is_straight PROPERTIES FIXTURES_REQUIRED band_history_runs)
# The broadband path carries no bands to write.
add_test(NAME band_history_needs_spectral COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -b tests/data/mixed_flush.csv -L warn --no-spectral --band_history integrated/band_history/broadband -o integrated/band_history_broadband.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(band_history_needs_spectral PROPERTIES WILL_FAIL TRUE)
# A point receiver, refracting through uniform air, ends its ray at the
# microphone as the straight path does: the two agree.
add_test(NAME refraction_uniform_point COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --ground_effect --atmosphere atmospheres/B407_uniform.csv --refraction linear --gridstep 492 -L warn -o integrated/refraction_uniform_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME straight_point COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --ground_effect --gridstep 492 -L warn -o integrated/straight_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(refraction_uniform_point straight_point PROPERTIES FIXTURES_SETUP refraction_point_runs)
add_test(NAME refraction_uniform_point_is_straight COMMAND golden_compare integrated/refraction_uniform_point.csv integrated/straight_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(refraction_uniform_point_is_straight PROPERTIES FIXTURES_REQUIRED refraction_point_runs)
# Real refraction, with wind, a shadow zone upwind, and each ground receiver,
# through the peak search against a march (-D).  The peak search bounds a
# refracted state by the same DirectionBounds as a straight one, taken along
# its launch direction at its ray tube's range; a skipped state must change
# nothing.  (The rays themselves are checked against traced ones in
# tests/test_refraction.cpp; eval_spectrum's asserts catch any non-finite level.)
add_test(NAME refraction_windy_point COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction linear --frame_bearing 210 --ground_effect --gridstep 492 -L warn -o integrated/refraction_windy_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME refraction_windy_point_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction linear --frame_bearing 210 --ground_effect --gridstep 492 -L warn -D integrated/refraction_windy_point_dnl.csv -o integrated/refraction_windy_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME refraction_windy_plate COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction linear --frame_bearing 210 --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --gridstep 492 -L warn -o integrated/refraction_windy_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME refraction_windy_plate_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction linear --frame_bearing 210 --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --gridstep 492 -L warn -D integrated/refraction_windy_plate_dnl.csv -o integrated/refraction_windy_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(refraction_windy_point refraction_windy_point_marched refraction_windy_plate refraction_windy_plate_marched PROPERTIES FIXTURES_SETUP refraction_windy_runs)
add_test(NAME stratified_windy_point COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction stratified --frame_bearing 210 --ground_effect --gridstep 984 -L warn -o integrated/stratified_windy_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_windy_point_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction stratified --frame_bearing 210 --ground_effect --gridstep 984 -L warn -D integrated/stratified_windy_point_dnl.csv -o integrated/stratified_windy_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_uniform_point COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --ground_effect --atmosphere atmospheres/B407_uniform.csv --refraction stratified --gridstep 492 -L warn -o integrated/stratified_uniform_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_windy_point stratified_windy_point_marched stratified_uniform_point PROPERTIES FIXTURES_SETUP stratified_runs)
add_test(NAME stratified_windowed_is_exact COMMAND golden_compare integrated/stratified_windy_point.csv integrated/stratified_windy_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_windowed_is_exact PROPERTIES FIXTURES_REQUIRED stratified_runs)
add_test(NAME stratified_uniform_point_is_straight COMMAND golden_compare integrated/stratified_uniform_point.csv integrated/straight_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_uniform_point_is_straight PROPERTIES FIXTURES_REQUIRED "stratified_runs;refraction_point_runs")
# A calm morning inversion capped near 350 ft (day 286's, run 286276, without its
# wind): a surface duct.  Past the reach of every direct ray only rays that bounced
# arrive (layered_arrivals), several per pair, and the peak search bounds their
# sum by the loudest one's direction over the harmonic sum of their ranges; a
# skipped state must change nothing.
add_test(NAME stratified_duct_point COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_inversion.csv --refraction stratified --ground_effect --gridstep 984 -L warn -o integrated/stratified_duct_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_duct_point_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_inversion.csv --refraction stratified --ground_effect --gridstep 984 -L warn -D integrated/stratified_duct_point_dnl.csv -o integrated/stratified_duct_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_duct_plate COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_inversion.csv --refraction stratified --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --gridstep 984 -L warn -o integrated/stratified_duct_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_duct_plate_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_inversion.csv --refraction stratified --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --gridstep 984 -L warn -D integrated/stratified_duct_plate_dnl.csv -o integrated/stratified_duct_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_duct_point stratified_duct_point_marched stratified_duct_plate stratified_duct_plate_marched PROPERTIES FIXTURES_SETUP stratified_duct_runs)
add_test(NAME stratified_duct_point_windowed_is_exact COMMAND golden_compare integrated/stratified_duct_point.csv integrated/stratified_duct_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_duct_plate_windowed_is_exact COMMAND golden_compare integrated/stratified_duct_plate.csv integrated/stratified_duct_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_duct_point_windowed_is_exact stratified_duct_plate_windowed_is_exact PROPERTIES FIXTURES_REQUIRED stratified_duct_runs)
# Pinned footprints, since windowed-versus-marched compares two runs of one binary and a change
# to the arrivals' summation, bounce loss, ground geometry or completeness passes it: dropping the
# extra arrivals from the energy sum and the primary's bounce loss moved the duct footprint by up
# to 7.8 dB with every stratified test green.  Regenerate with the build that changed them, after
# checking the change is meant (CLAUDE.md, "Every arrival, summed in energy").  Absolute 1e-4 dB:
# the root solves are not bit-stable across compilers.
add_test(NAME stratified_windy_point_golden COMMAND golden_compare integrated/stratified_windy_point.csv golden/stratified_windy_point.csv 1e-4 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_duct_point_golden COMMAND golden_compare integrated/stratified_duct_point.csv golden/stratified_duct_point.csv 1e-4 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_duct_plate_golden COMMAND golden_compare integrated/stratified_duct_plate.csv golden/stratified_duct_plate.csv 1e-4 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_windy_point_golden PROPERTIES FIXTURES_REQUIRED stratified_runs)
set_tests_properties(stratified_duct_point_golden stratified_duct_plate_golden PROPERTIES FIXTURES_REQUIRED stratified_duct_runs)
# --scattering fills shadows with scattered energy, a gain the peak search's bound has to allow: the
# creeping-wave and scattering caps were the survivors of mutation testing, with no CTest running the option.
add_test(NAME scattering_point COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction linear --frame_bearing 210 --ground_effect --scattering 3e-6 --gridstep 492 -L warn -o integrated/scattering_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME scattering_point_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction linear --frame_bearing 210 --ground_effect --scattering 3e-6 --gridstep 492 -L warn -D integrated/scattering_point_dnl.csv -o integrated/scattering_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME scattering_windowed_is_exact COMMAND golden_compare integrated/scattering_point.csv integrated/scattering_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# ... and the option does something: without it the footprint differs (golden_compare reports a mismatch, which is the pass; a missing file would say something else).
add_test(NAME scattering_is_not_off COMMAND golden_compare integrated/scattering_point.csv integrated/refraction_windy_point.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(scattering_is_not_off PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value" FIXTURES_REQUIRED "scattering_runs;refraction_windy_runs")
set_tests_properties(scattering_point scattering_point_marched PROPERTIES FIXTURES_SETUP scattering_runs)
set_tests_properties(scattering_windowed_is_exact PROPERTIES FIXTURES_REQUIRED scattering_runs)
add_test(NAME refraction_point_windowed_is_exact COMMAND golden_compare integrated/refraction_windy_point.csv integrated/refraction_windy_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME refraction_plate_windowed_is_exact COMMAND golden_compare integrated/refraction_windy_plate.csv integrated/refraction_windy_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(refraction_point_windowed_is_exact refraction_plate_windowed_is_exact PROPERTIES FIXTURES_REQUIRED refraction_windy_runs)
