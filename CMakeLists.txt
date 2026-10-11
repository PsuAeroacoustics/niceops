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
  src/pe_table.cpp
  src/airy.cpp
  src/exact_solutions.cpp
  src/weighting.cpp
  src/pnl.cpp
  src/band_frames.cpp
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
  src/condition_selection.cpp
  src/trajectory.cpp
  ${NICEOPS_PLOTTING_SOURCE}
  src/smoothing_spline.cpp
  src/manifest.cpp
  src/build_info.cpp
)
target_include_directories(niceops_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${NETCDF_INCLUDE_DIR}
)
# What the run manifest says built it.  The git state is read at every build, not at
# configure, which a later commit would leave stale; cmake/build_info.cmake rewrites the
# header only when the state changes, so an unchanged tree recompiles nothing.
set(NICEOPS_GENERATED_DIR ${CMAKE_CURRENT_BINARY_DIR}/generated)
file(MAKE_DIRECTORY ${NICEOPS_GENERATED_DIR})
add_custom_target(niceops_build_info
  COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
          -DOUTPUT=${NICEOPS_GENERATED_DIR}/niceops_build_info.hpp
          -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/build_info.cmake
  BYPRODUCTS ${NICEOPS_GENERATED_DIR}/niceops_build_info.hpp)
add_dependencies(niceops_core niceops_build_info)
target_include_directories(niceops_core PRIVATE ${NICEOPS_GENERATED_DIR})
set(NICEOPS_BUILD_INFO_DEFINITIONS
  "NICEOPS_COMPILER=\"${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}\""
  "NICEOPS_BUILD_TYPE=\"${CMAKE_BUILD_TYPE}\"")
if(NICEOPS_PLOTTING)
  list(APPEND NICEOPS_BUILD_INFO_DEFINITIONS NICEOPS_PLOTTING_ENABLED)
endif()
set_source_files_properties(src/build_info.cpp PROPERTIES COMPILE_DEFINITIONS "${NICEOPS_BUILD_INFO_DEFINITIONS}")
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
  tests/test_time_weighting.cpp
  tests/test_parallel.cpp
  tests/test_build.cpp
  tests/test_manifest.cpp
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
  tests/test_pe_table.cpp
  tests/test_benchmark.cpp
  tests/test_shadow_attenuation.cpp
  tests/test_shadow_transition.cpp
  tests/test_resolved_shadow.cpp
  tests/test_weighting.cpp
  tests/test_pnl.cpp
  tests/test_ground_effect.cpp
  tests/test_plate_receiver.cpp
  tests/test_spectral_propagation.cpp
  tests/test_predicates.cpp
  tests/test_delaunay.cpp
  tests/test_condition_selection.cpp
  tests/test_trajectory.cpp
  tests/test_csv.cpp
  tests/test_legacy_repair.cpp
  tests/test_smoothing_spline.cpp
)
target_link_libraries(niceops_unit_tests niceops_core Catch2::Catch2WithMain)

# Golden output comparison tool (numeric, elementwise)
add_executable(golden_compare tests/golden_compare.cpp)
# What --epnl_file and --band_frames wrote, checked against each other and the SEL.
add_executable(check_band_frames tests/check_band_frames.cpp)
target_link_libraries(check_band_frames niceops_core)
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
add_executable(check_tables_identical tests/check_tables_identical.cpp)
target_link_libraries(check_tables_identical niceops_core)
add_test(NAME pe_tools_thread_count_is_invisible COMMAND ${CMAKE_COMMAND}
  -DTOOL=$<TARGET_FILE:niceops_propagation_table>
  -DCOMPARE_TOOL=$<TARGET_FILE:niceops_pe_compare>
  -DCOMPARE=$<TARGET_FILE:check_tables_identical>
  -DOUT=${CMAKE_CURRENT_BINARY_DIR}/table_threads
  -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_table_threads.cmake
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_test(NAME propagation_table_reciprocity_and_turbulence COMMAND ${CMAKE_COMMAND}
  -DTOOL=$<TARGET_FILE:niceops_propagation_table>
  -DOUT=${CMAKE_CURRENT_BINARY_DIR}/table_turbulence
  -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_table_turbulence.cmake
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
# The load-time screen finds nothing in the golden database: the demo path never warns.
set_tests_properties(self_test PROPERTIES FAIL_REGULAR_EXPRESSION "dB from the median of (its|their)")
# It does find the S-76D fixed_load_factor file's 135-1450 dBA channels, names the worst,
# and with the mask drops them; on the golden database the mask changes nothing.
# PASS_REGULAR_EXPRESSION replaces the exit status, so the run is pinned by a plain test, the
# message by a second run, and the list by reading what the first left.
add_test(NAME screen_database_run COMMAND niceops -d databases/S-76D_fixed_load_factor.nod --gridstep 2000 -L warn -o integrated/screen_s76d.csv --screen_database integrated/screen_s76d_list.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(screen_database_run PROPERTIES FIXTURES_SETUP screen_s76d_run)
add_test(NAME screen_database_warns COMMAND niceops -d databases/S-76D_fixed_load_factor.nod --gridstep 2000 -L warn -o integrated/screen_s76d_warns.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(screen_database_warns PROPERTIES PASS_REGULAR_EXPRESSION "172 dBA cells in 4 of 164 spheres lie more than 40 dB from the median of their sphere.*1449\\.54")
add_test(NAME screen_database_lists_spikes COMMAND ${CMAKE_COMMAND} -E cat integrated/screen_s76d_list.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(screen_database_lists_spikes PROPERTIES FIXTURES_REQUIRED screen_s76d_run PASS_REGULAR_EXPRESSION "\n62,0\\.208[0-9]*,-7\\.77[0-9]*,0\\.00586[0-9]*,-120,180,,1449\\.5")
# On a copy: a guard that failed would replace the file.
add_test(NAME screen_database_victim COMMAND ${CMAKE_COMMAND} -E copy databases/AW139_fixed_load_factor.nod integrated/screen_victim.nod WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(screen_database_victim PROPERTIES FIXTURES_SETUP screen_victim)
add_test(NAME screen_database_refuses_overwriting_it COMMAND niceops -d integrated/screen_victim.nod --screen_database integrated/screen_victim.nod WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(screen_database_refuses_overwriting_it PROPERTIES FIXTURES_REQUIRED screen_victim PASS_REGULAR_EXPRESSION "would overwrite the database it reads")
add_test(NAME mask_suspect_cells_clean COMMAND niceops --no-spectral --no-ground_effect --mask_suspect_cells -L warn -o integrated/mask_suspect_clean.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(mask_suspect_cells_clean PROPERTIES FIXTURES_SETUP mask_clean_run)
add_test(NAME mask_suspect_cells_clean_is_golden COMMAND golden_compare integrated/mask_suspect_clean.csv golden/197421_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(mask_suspect_cells_clean_is_golden PROPERTIES FIXTURES_REQUIRED mask_clean_run)
add_test(NAME mask_suspect_cells_refuses_self_test COMMAND niceops --self_test --mask_suspect_cells WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(mask_suspect_cells_refuses_self_test PROPERTIES PASS_REGULAR_EXPRESSION "does not combine with --mask_suspect_cells")
add_test(NAME screen_database_refuses_two_databases COMMAND niceops -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197421AC.csv -d databases/S-76D_M3.nod -d databases/S-76D_M4.nod -o integrated/screen_two.csv --screen_database integrated/screen_two_list.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(screen_database_refuses_two_databases PROPERTIES PASS_REGULAR_EXPRESSION "run it once per database")
add_test(NAME unit_test COMMAND niceops_unit_tests WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# --manifest: the run's record.  The run is a plain test (exit status), the record is read back.
add_test(NAME manifest_run COMMAND niceops --no-spectral --no-ground_effect --gridstep 2000 -L warn -o integrated/manifest_demo.csv --manifest integrated/manifest_demo.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(manifest_run PROPERTIES FIXTURES_SETUP manifest_demo_run)
add_test(NAME manifest_records_run COMMAND ${CMAKE_COMMAND} -E cat integrated/manifest_demo.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(manifest_records_run PROPERTIES FIXTURES_REQUIRED manifest_demo_run PASS_REGULAR_EXPRESSION "\"role\": \"database\",[^}]*\"path\": \"databases/S-76D_M3\\.nod\",[^}]*\"sha256\": \"[0-9a-f]+\".*\"outputs\": \\[[^}]*\"path\": \"integrated/manifest_demo\\.csv\",[^}]*\"sha256\": \"[0-9a-f]+\".*\"exit_status\": 0")
# A run that stops still leaves its record, with the status it returned.
add_test(NAME manifest_records_failure COMMAND niceops --self_test --mask_suspect_cells --manifest integrated/manifest_failed.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(manifest_records_failure PROPERTIES FIXTURES_SETUP manifest_failed_run WILL_FAIL TRUE)
add_test(NAME manifest_records_failure_status COMMAND ${CMAKE_COMMAND} -E cat integrated/manifest_failed.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(manifest_records_failure_status PROPERTIES FIXTURES_REQUIRED manifest_failed_run PASS_REGULAR_EXPRESSION "\"errors\": 1,.*error: --self_test checks the database.*\"exit_status\": 2")
add_test(NAME manifest_victim COMMAND ${CMAKE_COMMAND} -E copy databases/S-76D.json integrated/manifest_victim.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(manifest_victim PROPERTIES FIXTURES_SETUP manifest_victim)
add_test(NAME manifest_refuses_overwriting_an_input COMMAND niceops -v integrated/manifest_victim.json --manifest integrated/manifest_victim.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(manifest_refuses_overwriting_an_input PROPERTIES FIXTURES_REQUIRED manifest_victim PASS_REGULAR_EXPRESSION "--manifest would overwrite the vehicle")
# --config: a file of options is the command line it spells.
add_test(NAME config_demo COMMAND niceops --config tests/data/demo_options.toml -L warn -o integrated/config_demo.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(config_demo PROPERTIES FIXTURES_SETUP config_demo_run)
add_test(NAME config_demo_is_golden COMMAND golden_compare integrated/config_demo.csv golden/197421_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(config_demo_is_golden PROPERTIES FIXTURES_REQUIRED config_demo_run)
add_test(NAME config_refuses_an_unknown_key COMMAND niceops --config tests/data/misspelled_options.toml -o integrated/config_misspelled.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(config_refuses_an_unknown_key PROPERTIES PASS_REGULAR_EXPRESSION "gridstpe")
# --provenance_file: a load-factor ceiling of 1.05 g holds 600 of the demo's 6551 states on B407.
# It needs every state, so the run marches, and the SEL must be the windowed run's.
add_test(NAME provenance_run COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 --gridstep 2000 --max_load_factor 1.05 -L error -o integrated/provenance_sel.csv --provenance_file integrated/provenance.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME provenance_windowed COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 --gridstep 2000 --max_load_factor 1.05 -L error -o integrated/provenance_windowed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(provenance_run provenance_windowed PROPERTIES FIXTURES_SETUP provenance_runs)
add_test(NAME provenance_leaves_sel_alone COMMAND golden_compare integrated/provenance_sel.csv integrated/provenance_windowed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME provenance_shares_held_states COMMAND ${CMAKE_COMMAND} -E cat integrated/provenance.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(provenance_leaves_sel_alone provenance_shares_held_states PROPERTIES FIXTURES_REQUIRED provenance_runs)
set_tests_properties(provenance_shares_held_states PROPERTIES PASS_REGULAR_EXPRESSION "exposure_db,share_mu_clamped,.*\n-4921\\.[0-9]+,-4921\\.[0-9]+,62\\.583[0-9]*,0,0,0,0\\.0955[0-9]*,0,31,0\\.485")
add_test(NAME provenance_refuses_ambient COMMAND niceops -a 40 --provenance_file integrated/provenance_ambient.csv -o integrated/provenance_ambient_sel.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(provenance_refuses_ambient PROPERTIES PASS_REGULAR_EXPRESSION "does not combine with -a")
add_test(NAME provenance_refuses_two_databases COMMAND niceops -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197421AC.csv -d databases/S-76D_M3.nod -d databases/S-76D_M4.nod -o integrated/provenance_two.csv --provenance_file integrated/provenance_two_list.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(provenance_refuses_two_databases PROPERTIES PASS_REGULAR_EXPRESSION "run it once per database")
add_test(NAME version_names_the_build COMMAND niceops --version WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(version_names_the_build PROPERTIES PASS_REGULAR_EXPRESSION "^[^ \n]+\n$")
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
# The period levels.  A noon flight lies wholly in the day: LAeq over the whole day is then Lden,
# whose day carries no penalty, and Lday is LAeq over 07:00-19:00.  A 01:00 flight lies in the
# night: Lnight is LAeq over 23:00-07:00, and Lday has no energy at all.
set(PERIOD_RUN -t tests/data/demo_geodetic_noclock.csv --origin 39.82914323 -74.40996141 --frame_bearing 210 --ground_elevation 54.2651 --gridstep 1640 -L warn)
add_test(NAME period_noon COMMAND niceops ${PERIOD_RUN} --start_utc 43200 --lden_file integrated/period_noon_lden.csv --lday_file integrated/period_noon_lday.csv --laeq_file integrated/period_noon_laeq24.csv --laeq_window 0 24 -o integrated/period_noon.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME period_noon_day_window COMMAND niceops ${PERIOD_RUN} --start_utc 43200 --laeq_file integrated/period_noon_laeq_day.csv --laeq_window 7 19 -o integrated/period_noon_b.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME period_night COMMAND niceops ${PERIOD_RUN} --start_utc 3600 --lnight_file integrated/period_night_lnight.csv --laeq_file integrated/period_night_laeq.csv --laeq_window 23 7 --lday_file integrated/period_night_lday.csv -o integrated/period_night.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(period_noon period_noon_day_window period_night PROPERTIES FIXTURES_SETUP period_runs)
add_test(NAME period_laeq_whole_day_is_lden_by_day COMMAND golden_compare integrated/period_noon_laeq24.csv integrated/period_noon_lden.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME period_lday_is_laeq_07_19 COMMAND golden_compare integrated/period_noon_lday.csv integrated/period_noon_laeq_day.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME period_lnight_is_laeq_23_07 COMMAND golden_compare integrated/period_night_lnight.csv integrated/period_night_laeq.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME period_lday_empty_at_night COMMAND ${CMAKE_COMMAND} -E cat integrated/period_night_lday.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(period_laeq_whole_day_is_lden_by_day period_lday_is_laeq_07_19 period_lnight_is_laeq_23_07 period_lday_empty_at_night PROPERTIES FIXTURES_REQUIRED period_runs)
# All seven rows of the 7 x 7 grid.
set_tests_properties(period_lday_empty_at_night PROPERTIES PASS_REGULAR_EXPRESSION "\n-inf,-inf,-inf,-inf,-inf,-inf,-inf\n-inf,-inf,-inf,-inf,-inf,-inf,-inf\n-inf,-inf,-inf,-inf,-inf,-inf,-inf\n-inf,-inf,-inf,-inf,-inf,-inf,-inf\n-inf,-inf,-inf,-inf,-inf,-inf,-inf\n-inf,-inf,-inf,-inf,-inf,-inf,-inf\n-inf,-inf,-inf,-inf,-inf,-inf,-inf\n")
add_test(NAME laeq_needs_a_window COMMAND niceops ${PERIOD_RUN} --start_utc 43200 --laeq_file integrated/laeq_no_window.csv -o integrated/laeq_no_window_sel.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(laeq_needs_a_window PROPERTIES PASS_REGULAR_EXPRESSION "--laeq_file needs --laeq_window")
add_test(NAME laeq_window_refuses_empty COMMAND niceops ${PERIOD_RUN} --start_utc 43200 --laeq_file integrated/laeq_empty.csv --laeq_window 7 7 -o integrated/laeq_empty_sel.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(laeq_window_refuses_empty PROPERTIES PASS_REGULAR_EXPRESSION "is empty; for a whole day give 0 24")
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
# Output metrics.  The demo track is flown 14:05-14:08 UTC, so --utc_offset puts
# it by day (0), at 22:30 local (8.4) and at 23:30 (9.4).  By day DNL, Lden and
# CNEL penalize nothing and agree; at 22:30 CNEL's night has begun (+10, as
# DNL's) and Lden's evening has not ended (+5); at 23:30 all three are night.
# Each pair below fails on a schedule whose boundary has moved.
set(METRICS_RUN -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --gridstep 1640 -L warn)
add_test(NAME metrics_by_day COMMAND niceops ${METRICS_RUN} -D integrated/metrics_day_dnl.csv --lden_file integrated/metrics_day_lden.csv --cnel_file integrated/metrics_day_cnel.csv -o integrated/metrics_day.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME metrics_at_22h30 COMMAND niceops ${METRICS_RUN} --utc_offset 8.4 -D integrated/metrics_2230_dnl.csv --lden_file integrated/metrics_2230_lden.csv --cnel_file integrated/metrics_2230_cnel.csv -o integrated/metrics_2230.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME metrics_at_23h30 COMMAND niceops ${METRICS_RUN} --utc_offset 9.4 -D integrated/metrics_2330_dnl.csv --lden_file integrated/metrics_2330_lden.csv --cnel_file integrated/metrics_2330_cnel.csv -o integrated/metrics_2330.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# --full_history governs DNL too: each operation's whole history rather than its
# own 10 dB-down window, so the two must differ.
add_test(NAME metrics_full_history COMMAND niceops ${METRICS_RUN} --full_history -D integrated/metrics_full_dnl.csv -o integrated/metrics_full.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# C-SEL: limited to the 1 kHz band, where the C and A curves are both 0 dB (to
# 1.5e-6 dB), it is SEL, which checks the window, the reception weights and the
# per-state bands; over every band it is not.
add_test(NAME metrics_csel_1khz COMMAND niceops ${METRICS_RUN} --band_limits 900 1100 --csel_file integrated/metrics_csel_1khz.csv -o integrated/metrics_sel_1khz.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME metrics_csel COMMAND niceops ${METRICS_RUN} --csel_file integrated/metrics_csel.csv -o integrated/metrics_sel.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# LAmax: the fast weighting is not the instantaneous peak.
add_test(NAME metrics_lamax_fast COMMAND niceops ${METRICS_RUN} --lamax_file integrated/metrics_lafmax.csv --nax_file integrated/metrics_nax.csv --nax_threshold 45 --ta_file integrated/metrics_ta.csv --ta_threshold 40 -o integrated/metrics_lamax_run.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME metrics_lamax_instantaneous COMMAND niceops ${METRICS_RUN} --lamax_weighting instantaneous --lamax_file integrated/metrics_lamax_peak.csv -o integrated/metrics_lamax_peak_run.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(metrics_by_day metrics_at_22h30 metrics_at_23h30 metrics_full_history metrics_csel_1khz metrics_csel metrics_lamax_fast metrics_lamax_instantaneous PROPERTIES FIXTURES_SETUP metrics_runs)
add_test(NAME lden_by_day_is_dnl COMMAND golden_compare integrated/metrics_day_lden.csv integrated/metrics_day_dnl.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME cnel_by_day_is_dnl COMMAND golden_compare integrated/metrics_day_cnel.csv integrated/metrics_day_dnl.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME cnel_at_22h30_is_dnl COMMAND golden_compare integrated/metrics_2230_cnel.csv integrated/metrics_2230_dnl.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME lden_at_23h30_is_dnl COMMAND golden_compare integrated/metrics_2330_lden.csv integrated/metrics_2330_dnl.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME csel_at_1khz_is_sel COMMAND golden_compare integrated/metrics_csel_1khz.csv integrated/metrics_sel_1khz.csv 1e-5 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# The "differs" cases pass on golden_compare's mismatch message, not on its exit
# status, so a missing file cannot satisfy them.
add_test(NAME lden_at_22h30_is_not_dnl COMMAND golden_compare integrated/metrics_2230_lden.csv integrated/metrics_2230_dnl.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME dnl_full_history_is_not_windowed COMMAND golden_compare integrated/metrics_full_dnl.csv integrated/metrics_day_dnl.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME csel_is_not_sel COMMAND golden_compare integrated/metrics_csel.csv integrated/metrics_sel.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME lafmax_is_not_the_peak COMMAND golden_compare integrated/metrics_lafmax.csv integrated/metrics_lamax_peak.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(lden_by_day_is_dnl cnel_by_day_is_dnl cnel_at_22h30_is_dnl lden_at_23h30_is_dnl csel_at_1khz_is_sel lden_at_22h30_is_not_dnl dnl_full_history_is_not_windowed csel_is_not_sel lafmax_is_not_the_peak PROPERTIES FIXTURES_REQUIRED metrics_runs)
set_tests_properties(lden_at_22h30_is_not_dnl dnl_full_history_is_not_windowed csel_is_not_sel lafmax_is_not_the_peak PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value")
# The call site: LAFmax, NA64 and TA70 on the demo case's broadband database (the
# one regression_golden pins), against goldens checked independently -- a brute-
# force fast detector on a 5e-5 s reception-time grid built from the -S state and
# the -i histories agreed to 1e-4 dB, where the same in emission time is 0.018 dB
# off.  So these fail if the weighting hears emission time, or NAx counts the
# instantaneous peak (64 dBA sits between the two at three cells).
add_test(NAME metrics_broadband COMMAND niceops --gridstep 1640 -L warn --lamax_file integrated/metrics_bb_lafmax.csv --nax_file integrated/metrics_bb_na64.csv --nax_threshold 64 --ta_file integrated/metrics_bb_ta70.csv --ta_threshold 70 -o integrated/metrics_bb.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(metrics_broadband PROPERTIES FIXTURES_SETUP metrics_broadband_run)
add_test(NAME lafmax_golden COMMAND golden_compare integrated/metrics_bb_lafmax.csv golden/metrics_197421_lafmax.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME nax_golden COMMAND golden_compare integrated/metrics_bb_na64.csv golden/metrics_197421_na64.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ta_golden COMMAND golden_compare integrated/metrics_bb_ta70.csv golden/metrics_197421_ta70.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(lafmax_golden nax_golden ta_golden PROPERTIES FIXTURES_REQUIRED metrics_broadband_run)
# What the metrics refuse.  A 1 s geodetic track is resampled onto the 0.5 s -s
# default, too coarse for the fast detector; at -s 0.05 it runs.
set(METRICS_GEODETIC -t tests/data/demo_geodetic_noclock.csv --origin 39.82914323 -74.40996141 --frame_bearing 210 --ground_elevation 54.2651 --gridstep 1640)
add_test(NAME lamax_weighting_refuses_a_coarse_history COMMAND niceops ${METRICS_GEODETIC} --lamax_file integrated/coarse_lamax.csv -o integrated/coarse_lamax_run.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(lamax_weighting_refuses_a_coarse_history PROPERTIES PASS_REGULAR_EXPRESSION "needs the level history sampled at least every 0.0625 s")
add_test(NAME lamax_weighting_runs_on_a_finer_step COMMAND niceops ${METRICS_GEODETIC} -s 0.05 --lamax_file integrated/fine_lamax.csv -o integrated/fine_lamax_run.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(lamax_weighting_runs_on_a_finer_step PROPERTIES FAIL_REGULAR_EXPRESSION "needs the level history")
add_test(NAME lden_needs_a_clock COMMAND niceops ${METRICS_GEODETIC} --lden_file integrated/noclock_lden.csv -o integrated/noclock_lden_run.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(lden_needs_a_clock PROPERTIES PASS_REGULAR_EXPRESSION "--lden_file needs the time of day")
add_test(NAME lamax_weighting_needs_an_output COMMAND niceops ${METRICS_RUN} --lamax_weighting slow -o integrated/weighting_alone.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(lamax_weighting_needs_an_output PROPERTIES PASS_REGULAR_EXPRESSION "none of them is given")
add_test(NAME ta_needs_a_threshold COMMAND niceops ${METRICS_RUN} --ta_file integrated/ta_unset.csv -o integrated/ta_unset_run.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ta_needs_a_threshold PROPERTIES PASS_REGULAR_EXPRESSION "--ta_file needs --ta_threshold")
add_test(NAME csel_needs_spectra COMMAND niceops --gridstep 1640 --csel_file integrated/csel_broadband.csv -o integrated/csel_broadband_run.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(csel_needs_spectra PROPERTIES PASS_REGULAR_EXPRESSION "needs a database carrying per-band spectra, to C-weight")
add_test(NAME csel_refuses_no_spectral COMMAND niceops ${METRICS_RUN} --no-spectral --csel_file integrated/csel_refused.csv -o integrated/csel_refused_run.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(csel_refuses_no_spectral PROPERTIES PASS_REGULAR_EXPRESSION "--csel_file needs the spectral path")
add_test(NAME csel_refuses_db_cutoff COMMAND niceops ${METRICS_RUN} --db_cutoff 40 --csel_file integrated/csel_refused.csv -o integrated/csel_refused_run.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(csel_refuses_db_cutoff PROPERTIES PASS_REGULAR_EXPRESSION "does not combine with --db_cutoff")
add_test(NAME csel_refuses_ambient COMMAND niceops ${METRICS_RUN} -a 30 --csel_file integrated/csel_refused.csv -o integrated/csel_refused_run.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(csel_refuses_ambient PROPERTIES PASS_REGULAR_EXPRESSION "does not combine with -a")
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
# Plate microphones read the table at its lowest receiver height plus K, the plate over bare ground
# there at the pair's elevation and azimuth.  Through the ground-model table that must be the plate
# run (0.13 dB at worst on this grid, through its 1.2 m row; K taken at 2 degrees for every pair
# would have been 5.6 dB off), and K enters the peak search's bound, so windowed equals marched.
# A table without a row at 0 is warned about.
add_test(NAME propagation_table_plate COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --flow_resistance 225 --ground_effect --plate_table plates/gr1425_db225.nc --propagation_table tables/ground_model_test.nc --gridstep 492 -L warn -o integrated/propagation_table_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_plate_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --flow_resistance 225 --ground_effect --plate_table plates/gr1425_db225.nc --propagation_table tables/ground_model_test.nc --gridstep 492 -L warn -D integrated/propagation_table_plate_dnl.csv -o integrated/propagation_table_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_plate_reference COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --flow_resistance 225 --ground_effect --exact_ground --plate_table plates/gr1425_db225.nc --gridstep 492 -L warn -o integrated/propagation_table_plate_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_plate propagation_table_plate_marched PROPERTIES FIXTURES_REQUIRED propagation_table_file FIXTURES_SETUP propagation_table_plate_runs)
set_tests_properties(propagation_table_plate_reference PROPERTIES FIXTURES_SETUP propagation_table_plate_runs)
add_test(NAME propagation_table_plate_windowed_is_exact COMMAND golden_compare integrated/propagation_table_plate.csv integrated/propagation_table_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_plate_is_the_plate_model COMMAND golden_compare integrated/propagation_table_plate.csv integrated/propagation_table_plate_reference.csv 0.15 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_plate_windowed_is_exact propagation_table_plate_is_the_plate_model PROPERTIES FIXTURES_REQUIRED propagation_table_plate_runs)
add_test(NAME propagation_table_plate_warns_raised_row COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --flow_resistance 225 --ground_effect --plate_table plates/gr1425_db225.nc --propagation_table tables/ground_model_test.nc --gridstep 4000 -L warn -o integrated/propagation_table_plate_warns.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_plate_warns_raised_row PROPERTIES FIXTURES_REQUIRED propagation_table_file PASS_REGULAR_EXPRESSION "at its lowest receiver height, 3\\.937.*Give the table a receiver height of 0")
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
# --blend and --load_factor_correction (condition_selection.hpp).  Barycentric is the default, so
# naming it, and the default one-sided law, must reproduce the golden; nearest and aam must run and
# move it (golden_compare's mismatch is the pass), and differ from each other.  The demo database is
# S-76D_M3, whose thrust axis is spheres, so aam there searches its thrust layers.
add_test(NAME blend_barycentric COMMAND niceops --no-spectral --no-ground_effect --blend barycentric --load_factor_correction one_sided -L warn -o integrated/blend_barycentric.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME blend_nearest COMMAND niceops --no-spectral --no-ground_effect --blend nearest -L warn -o integrated/blend_nearest.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME blend_aam COMMAND niceops --no-spectral --no-ground_effect --blend aam -L warn -o integrated/blend_aam.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(blend_barycentric blend_nearest blend_aam PROPERTIES FIXTURES_SETUP blend_runs)
add_test(NAME blend_barycentric_is_golden COMMAND golden_compare integrated/blend_barycentric.csv golden/197421_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME blend_nearest_moves_footprint COMMAND golden_compare integrated/blend_nearest.csv golden/197421_LAE.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME blend_aam_moves_footprint COMMAND golden_compare integrated/blend_aam.csv golden/197421_LAE.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME blend_aam_is_not_nearest COMMAND golden_compare integrated/blend_aam.csv integrated/blend_nearest.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(blend_nearest_moves_footprint blend_aam_moves_footprint blend_aam_is_not_nearest PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value")
set_tests_properties(blend_barycentric_is_golden blend_nearest_moves_footprint blend_aam_moves_footprint blend_aam_is_not_nearest PROPERTIES FIXTURES_REQUIRED blend_runs)
# The spectral path's windowed peak search bounds every state by its blend's spheres and its
# load-factor shift, so each policy, and the two-sided law's negative shifts, must leave it exact
# against a march (-D forces the march) on a fixed_load_factor database.
add_test(NAME blend_aam_spectral COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --blend aam --load_factor_correction two_sided --gridstep 492 -L warn -o integrated/blend_aam_spectral.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME blend_aam_spectral_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --blend aam --load_factor_correction two_sided --gridstep 492 -L warn -D integrated/blend_aam_spectral_dnl.csv -o integrated/blend_aam_spectral_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME blend_nearest_spectral COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --blend nearest --gridstep 492 -L warn -o integrated/blend_nearest_spectral.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME blend_nearest_spectral_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --blend nearest --gridstep 492 -L warn -D integrated/blend_nearest_spectral_dnl.csv -o integrated/blend_nearest_spectral_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(blend_aam_spectral blend_aam_spectral_marched blend_nearest_spectral blend_nearest_spectral_marched PROPERTIES FIXTURES_SETUP blend_spectral_runs)
add_test(NAME blend_aam_windowed_is_exact COMMAND golden_compare integrated/blend_aam_spectral.csv integrated/blend_aam_spectral_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME blend_nearest_windowed_is_exact COMMAND golden_compare integrated/blend_nearest_spectral.csv integrated/blend_nearest_spectral_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# ... and the two-sided law does something on this track, whose load factor runs both sides of 1 g:
# atmosphere_homogeneous is the same run, barycentric and one-sided.
add_test(NAME load_factor_two_sided COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --load_factor_correction two_sided --min_load_factor 0.5 --gridstep 492 -L warn -o integrated/load_factor_two_sided.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(load_factor_two_sided PROPERTIES FIXTURES_SETUP blend_spectral_runs)
add_test(NAME load_factor_two_sided_moves_footprint COMMAND golden_compare integrated/load_factor_two_sided.csv integrated/atmosphere_homogeneous.csv 0.001 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(load_factor_two_sided_moves_footprint PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value")
set_tests_properties(blend_aam_windowed_is_exact blend_nearest_windowed_is_exact load_factor_two_sided_moves_footprint PROPERTIES FIXTURES_REQUIRED "blend_spectral_runs;atmosphere_runs")
# Refused, on the message: the floor without the two-sided law, and the law on a database whose load
# factors are spheres.
add_test(NAME min_load_factor_needs_two_sided COMMAND niceops --min_load_factor 0.3 -L warn -o integrated/min_load_factor_needs_two_sided.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(min_load_factor_needs_two_sided PROPERTIES PASS_REGULAR_EXPRESSION "needs --load_factor_correction two_sided")
add_test(NAME two_sided_needs_fixed_load_factor COMMAND niceops --no-spectral --no-ground_effect --load_factor_correction two_sided --gridstep 2000 -L warn -o integrated/two_sided_needs_fixed_load_factor.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(two_sided_needs_fixed_load_factor PROPERTIES PASS_REGULAR_EXPRESSION "is not one: its load factors are spheres")
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
# --epnl_file and --band_frames (src/pnl.hpp, src/band_frames.hpp).  The run writes both over
# flush plate microphones and pole microphones; check_band_frames recomputes every EPNL from the
# netCDF frames, checks the files against each other and for plausibility, and holds the frames'
# A-weighted energy to the --full_history SEL of the same case: every state's energy, in reception
# time, is in the frames.  The frames force a march, which must leave the footprint equal to the
# peak search's.  A coarse grid checks the grid layout.
add_test(NAME epnl_run COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -b tests/data/mixed_flush.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc -L warn --epnl_file integrated/epnl/mixed_EPNL.csv --band_frames integrated/epnl/mixed_frames.nc -o integrated/epnl/mixed_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME epnl_full_history COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -b tests/data/mixed_flush.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc -L warn --full_history -o integrated/epnl/mixed_full_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME epnl_windowed COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -b tests/data/mixed_flush.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc -L warn -o integrated/epnl/mixed_windowed_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME epnl_grid COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --gridstep 1312 -L warn --epnl_file integrated/epnl/grid_EPNL.csv --band_frames integrated/epnl/grid_frames.nc -o integrated/epnl/grid_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# The same case's band histories, from which check_band_frames --history rebuilds every frame: the
# energy check cannot see where in time the energy went, and a state filed at its emission time
# (or at t - r/c) passed it.  And the case once more with --no-epnl_masking, the tone correction as
# written: on these microphones inaudible bands give tones (one of the eight PNLTM records, and
# others inside the durations), so some EPNL must read higher.
add_test(NAME epnl_band_history COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -b tests/data/mixed_flush.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc -L warn --band_history integrated/epnl/history/mixed -o integrated/epnl/mixed_history_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME epnl_unmasked COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -b tests/data/mixed_flush.csv --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc -L warn --no-epnl_masking --epnl_file integrated/epnl/mixed_unmasked_EPNL.csv --band_frames integrated/epnl/mixed_unmasked_frames.nc -o integrated/epnl/mixed_unmasked_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(epnl_run epnl_full_history epnl_windowed epnl_grid epnl_band_history epnl_unmasked PROPERTIES FIXTURES_SETUP epnl_runs)
add_test(NAME epnl_frames_follow_reception_time COMMAND check_band_frames --history integrated/epnl/history/mixed integrated/epnl/mixed_frames.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME epnl_unmasked_is_consistent COMMAND check_band_frames --no-masking integrated/epnl/mixed_unmasked_frames.nc integrated/epnl/mixed_unmasked_EPNL.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME epnl_masking_lowers_tones COMMAND check_band_frames --compare integrated/epnl/mixed_EPNL.csv integrated/epnl/mixed_unmasked_EPNL.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(epnl_frames_follow_reception_time epnl_unmasked_is_consistent epnl_masking_lowers_tones PROPERTIES FIXTURES_REQUIRED epnl_runs)
add_test(NAME epnl_frames_are_consistent COMMAND check_band_frames integrated/epnl/mixed_frames.nc integrated/epnl/mixed_EPNL.csv integrated/epnl/mixed_full_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME epnl_grid_frames_are_consistent COMMAND check_band_frames integrated/epnl/grid_frames.nc integrated/epnl/grid_EPNL.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME epnl_leaves_footprint COMMAND golden_compare integrated/epnl/mixed_LAE.csv integrated/epnl/mixed_windowed_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(epnl_frames_are_consistent epnl_grid_frames_are_consistent epnl_leaves_footprint PROPERTIES FIXTURES_REQUIRED epnl_runs)
# Refused, each with its reason: no bands on the broadband path or in a database without spectra,
# an ambient with no spectrum, and band limits that leave out a PNL band.
add_test(NAME epnl_needs_spectral COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -b tests/data/mixed_flush.csv -L warn --no-spectral --epnl_file integrated/epnl/refused.csv -o integrated/epnl/refused_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(epnl_needs_spectral PROPERTIES PASS_REGULAR_EXPRESSION "--epnl_file needs the spectral path")
add_test(NAME band_frames_need_spectra COMMAND niceops -d databases/S-76D_M3.nod -t tracks/S-76D/197421AC.csv -b tests/data/mixed_flush.csv -L warn --band_frames integrated/epnl/refused.nc -o integrated/epnl/refused_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(band_frames_need_spectra PROPERTIES PASS_REGULAR_EXPRESSION "--band_frames needs a database carrying per-band spectra")
add_test(NAME epnl_refuses_ambient COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -b tests/data/mixed_flush.csv -L warn -a 40 --epnl_file integrated/epnl/refused.csv -o integrated/epnl/refused_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(epnl_refuses_ambient PROPERTIES PASS_REGULAR_EXPRESSION "--epnl_file does not combine with -a")
add_test(NAME epnl_needs_its_bands COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -b tests/data/mixed_flush.csv -L warn --band_limits 100 20000 --epnl_file integrated/epnl/refused.csv -o integrated/epnl/refused_LAE.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(epnl_needs_its_bands PROPERTIES PASS_REGULAR_EXPRESSION "leaves out the 50 Hz band")
# Several trajectories and one name: one file per track, named after it.
add_test(NAME epnl_multi_trajectory_run COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197405AC.csv -b tests/data/mixed_flush.csv -L warn --epnl_file integrated/epnl_multi/e.csv --band_frames integrated/epnl_multi/f.nc -o integrated/epnl_multi/lae.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME epnl_multi_trajectory_names COMMAND ${CMAKE_COMMAND} "-DFILES=integrated/epnl_multi/197421AC_EPNL.csv;integrated/epnl_multi/197405AC_EPNL.csv;integrated/epnl_multi/197421AC_band_frames.nc;integrated/epnl_multi/197405AC_band_frames.nc" "-DABSENT=integrated/epnl_multi/e.csv;integrated/epnl_multi/f.nc" -P tests/check_files.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(epnl_multi_trajectory_run PROPERTIES FIXTURES_SETUP epnl_multi_runs)
set_tests_properties(epnl_multi_trajectory_names PROPERTIES FIXTURES_REQUIRED epnl_multi_runs)
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
# --shadow_model resolved: each band's shadow boundary from the near-ground profile it resolves, and
# cramped bands' field from the residue series at the source's own height (or the straight path's
# ground term).  The peak search must stay exact with it on: the series and the table are capped at
# the transition's largest power ratio, and the straight ground term at the ground's own bound.
# Twins of stratified_windowed_is_exact (point and plate), of the duct's, and of the straight mode.
add_test(NAME stratified_resolved_point COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction stratified --frame_bearing 210 --ground_effect --shadow_model resolved --gridstep 984 -L warn -o integrated/stratified_resolved_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_resolved_point_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction stratified --frame_bearing 210 --ground_effect --shadow_model resolved --gridstep 984 -L warn -D integrated/stratified_resolved_point_dnl.csv -o integrated/stratified_resolved_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_resolved_plate COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction stratified --frame_bearing 210 --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --shadow_model resolved --gridstep 984 -L warn -o integrated/stratified_resolved_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_resolved_plate_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction stratified --frame_bearing 210 --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --shadow_model resolved --gridstep 984 -L warn -D integrated/stratified_resolved_plate_dnl.csv -o integrated/stratified_resolved_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_resolved_straight_point COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction stratified --frame_bearing 210 --ground_effect --cramped_ground straight --shadow_model resolved --gridstep 984 -L warn -o integrated/stratified_resolved_straight_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_resolved_straight_point_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction stratified --frame_bearing 210 --ground_effect --cramped_ground straight --shadow_model resolved --gridstep 984 -L warn -D integrated/stratified_resolved_straight_point_dnl.csv -o integrated/stratified_resolved_straight_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_resolved_duct_point COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_inversion.csv --refraction stratified --ground_effect --shadow_model resolved --gridstep 984 -L warn -o integrated/stratified_resolved_duct_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_resolved_duct_point_marched COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_inversion.csv --refraction stratified --ground_effect --shadow_model resolved --gridstep 984 -L warn -D integrated/stratified_resolved_duct_point_dnl.csv -o integrated/stratified_resolved_duct_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_resolved_point stratified_resolved_point_marched stratified_resolved_plate stratified_resolved_plate_marched stratified_resolved_straight_point stratified_resolved_straight_point_marched stratified_resolved_duct_point stratified_resolved_duct_point_marched PROPERTIES FIXTURES_SETUP stratified_resolved_runs)
add_test(NAME stratified_resolved_point_windowed_is_exact COMMAND golden_compare integrated/stratified_resolved_point.csv integrated/stratified_resolved_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_resolved_plate_windowed_is_exact COMMAND golden_compare integrated/stratified_resolved_plate.csv integrated/stratified_resolved_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_resolved_straight_point_windowed_is_exact COMMAND golden_compare integrated/stratified_resolved_straight_point.csv integrated/stratified_resolved_straight_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_resolved_duct_point_windowed_is_exact COMMAND golden_compare integrated/stratified_resolved_duct_point.csv integrated/stratified_resolved_duct_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_resolved_point_windowed_is_exact stratified_resolved_plate_windowed_is_exact stratified_resolved_straight_point_windowed_is_exact stratified_resolved_duct_point_windowed_is_exact PROPERTIES FIXTURES_REQUIRED stratified_resolved_runs)
# ... and it does something (a mismatch is the pass), while the default is the legacy transition, untouched.
add_test(NAME stratified_resolved_is_not_legacy COMMAND golden_compare integrated/stratified_resolved_point.csv integrated/stratified_windy_point.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_resolved_is_not_legacy PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value" FIXTURES_REQUIRED "stratified_resolved_runs;stratified_runs")
add_test(NAME shadow_model_needs_stratified COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction linear --shadow_model resolved --gridstep 1640 -o integrated/shadow_model_needs_stratified.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(shadow_model_needs_stratified PROPERTIES PASS_REGULAR_EXPRESSION "--shadow_model resolved sets the shadow boundaries --refraction stratified casts")
add_test(NAME shadow_phase_needs_resolved COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction stratified --shadow_phase 0.5 --gridstep 1640 -o integrated/shadow_phase_needs_resolved.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(shadow_phase_needs_resolved PROPERTIES PASS_REGULAR_EXPRESSION "--shadow_phase and --cramped_ground shape --shadow_model resolved")
# Discretization independence: example_windy.csv with a calm, isothermal floor 0.1 ft thick (a level
# at 0.1 ft, as fast as the ground) against the same profile without it.  Ray theory sent the upwind
# shadow's limiting ray along the floor without limit, and flush microphones on it came out up to
# 39 dB low on this grid; no band of the run resolves a floor that thin (resolved_pass, stratified.cpp).
# What may remain is the floor itself, which reshapes the lowest layer (0.05 dB at 1.2 m), and the
# cap's slack over what ray theory travels in a resolved turning, (1.5 - 4/3) k h^2, which leaves the
# boundary under 1% farther out at 10 kHz: 0.12 dB at worst on the plates, at the shadow's edge.
# 0.3 dB absolute; the unfixed solve misses by 39 dB.
add_test(NAME stratified_windy_plate COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction stratified --frame_bearing 210 --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --gridstep 984 -L warn -o integrated/stratified_windy_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_calm_floor_plate COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy_calm_floor.csv --refraction stratified --frame_bearing 210 --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --gridstep 984 -L warn -o integrated/stratified_calm_floor_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_calm_floor_point COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy_calm_floor.csv --refraction stratified --frame_bearing 210 --ground_effect --gridstep 984 -L warn -o integrated/stratified_calm_floor_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_windy_plate stratified_calm_floor_plate stratified_calm_floor_point PROPERTIES FIXTURES_SETUP stratified_calm_floor_runs)
add_test(NAME stratified_calm_floor_plate_is_floorless COMMAND golden_compare integrated/stratified_calm_floor_plate.csv integrated/stratified_windy_plate.csv 0.3 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_calm_floor_point_is_floorless COMMAND golden_compare integrated/stratified_calm_floor_point.csv integrated/stratified_windy_point.csv 0.3 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_calm_floor_plate_is_floorless PROPERTIES FIXTURES_REQUIRED stratified_calm_floor_runs)
set_tests_properties(stratified_calm_floor_point_is_floorless PROPERTIES FIXTURES_REQUIRED "stratified_calm_floor_runs;stratified_runs")
# The run's highest band reaches the layered solve (process_trajectory hands LayeredSoundSpeed its top
# band, path_layers each path's Layers).  With --band_limits ending at the 200 Hz band the floor is
# followed farther than at the 10 kHz default, and dropping the band from that call moves these plates;
# nothing else runs the solve at another frequency.  Pinned at 1e-4 dB like the
# other stratified goldens; regenerate with the build that changes the cap on purpose.
add_test(NAME stratified_calm_floor_plate_200hz COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy_calm_floor.csv --refraction stratified --frame_bearing 210 --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc --gridstep 984 --band_limits 20 224 -L warn -o integrated/stratified_calm_floor_plate_200hz.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_calm_floor_plate_200hz PROPERTIES FIXTURES_SETUP stratified_calm_floor_200hz_run)
add_test(NAME stratified_calm_floor_plate_200hz_golden COMMAND golden_compare integrated/stratified_calm_floor_plate_200hz.csv golden/stratified_calm_floor_plate_200hz.csv 1e-4 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_calm_floor_plate_200hz_golden PROPERTIES FIXTURES_REQUIRED stratified_calm_floor_200hz_run)
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
# --wind builds the profile in memory from one surface wind.  The run has to be exactly the one fed
# that profile through --atmosphere (--write_atmosphere writes it), and refining its levels from 16
# a decade to 64 must not move a plate microphone by 0.1 dB (it moves them 0.016).  The plate
# over the Coyle array under stratified refraction is the configuration measured against the 2017
# data, and the one whose thin-layer defect a calm level above the ground triggers.
set(WIND_RUN -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv -b arrays/coyle_approach.csv --refraction stratified --frame_bearing 210 --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc -L warn)
add_test(NAME wind_built COMMAND niceops ${WIND_RUN} --humidity 30 --wind 10kt,270,33 --write_atmosphere integrated/wind_profile.csv -o integrated/wind_built.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME wind_refined COMMAND niceops ${WIND_RUN} --humidity 30 --wind 10kt,270,33 --profile_levels 64 -o integrated/wind_refined.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(wind_built wind_refined PROPERTIES FIXTURES_SETUP wind_runs)
add_test(NAME wind_replayed COMMAND niceops ${WIND_RUN} --atmosphere integrated/wind_profile.csv -o integrated/wind_replayed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(wind_replayed PROPERTIES FIXTURES_REQUIRED wind_runs FIXTURES_SETUP wind_replay_runs)
add_test(NAME wind_is_its_written_profile COMMAND golden_compare integrated/wind_replayed.csv integrated/wind_built.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(wind_is_its_written_profile PROPERTIES FIXTURES_REQUIRED wind_replay_runs)
add_test(NAME wind_is_discretization_independent COMMAND golden_compare integrated/wind_refined.csv integrated/wind_built.csv 0.1 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(wind_is_discretization_independent PROPERTIES FIXTURES_REQUIRED wind_runs)
# A calm observation is still air, stated rather than assumed, so a database referred to the air
# takes the ground velocity as the airspeed instead of refusing for want of wind.
add_test(NAME wind_calm COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 30 -t tracks/S-76D/197421AC.csv --wind 0kt,0 --gridstep 1640 -L warn -o integrated/wind_calm.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME wind_calm_air COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 30 -t tracks/S-76D/197421AC.csv --wind 0kt,0 --speed_reference air --gridstep 1640 -L warn -o integrated/wind_calm_air.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(wind_calm wind_calm_air PROPERTIES FIXTURES_SETUP wind_calm_runs)
add_test(NAME wind_calm_air_is_ground COMMAND golden_compare integrated/wind_calm_air.csv integrated/wind_calm.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(wind_calm_air_is_ground PROPERTIES FIXTURES_REQUIRED wind_calm_runs)
# The written profile against its closed forms, through the command line: the log law, the dry
# lapse, the barometric formula for it and the held water vapor, computed independently
# (tests/data/wind_profile_expected.py) for a run in feet that sets every shaping option and all
# three surface overrides.  The replay test above compares the run with itself and cannot see a
# length left in feet, an override not applied or --profile_levels not passed on; this can.
add_test(NAME wind_profile_written COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --gridstep 1640 -L warn --wind 12kt,250,6.5 --roughness_length 0.1 --surface_layer_depth 3000 --profile_levels 10 --temperature 288.15 --humidity 40 --pressure 100 --write_atmosphere integrated/wind_profile_law.csv -o integrated/wind_profile_law_levels.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(wind_profile_written PROPERTIES FIXTURES_SETUP wind_law_runs)
add_test(NAME wind_profile_is_its_law COMMAND golden_compare integrated/wind_profile_law.csv tests/data/wind_profile_expected.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(wind_profile_is_its_law PROPERTIES FIXTURES_REQUIRED wind_law_runs)
# What --wind refuses, by message: a speed without its unit, a file profile beside it, its shaping
# options without it, and a surface humidity nobody gave.  The unvalidated classes run, and say so.
add_test(NAME wind_speed_needs_unit COMMAND niceops --wind 10,270 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(wind_speed_needs_unit PROPERTIES PASS_REGULAR_EXPRESSION "--wind: name the speed's unit")
add_test(NAME wind_with_atmosphere_refused COMMAND niceops --wind 10kt,270 --atmosphere atmospheres/isa.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(wind_with_atmosphere_refused PROPERTIES PASS_REGULAR_EXPRESSION "give one or the other")
add_test(NAME wind_options_need_wind COMMAND niceops --roughness_length 0.1 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(wind_options_need_wind PROPERTIES PASS_REGULAR_EXPRESSION "--roughness_length shapes the profile --wind builds, and there is no --wind")
add_test(NAME wind_anemometer_below_roughness_refused COMMAND niceops --wind 10kt,270,0.05 --roughness_length 0.1 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(wind_anemometer_below_roughness_refused PROPERTIES PASS_REGULAR_EXPRESSION "the anemometer \\(0.05 ft\\) has to stand above it")
add_test(NAME wind_needs_humidity COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --wind 10kt,270 --gridstep 1640 -o integrated/wind_needs_humidity.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(wind_needs_humidity PROPERTIES PASS_REGULAR_EXPRESSION "--wind builds the humidity profile from the surface humidity, and none is set")
add_test(NAME wind_stability_class_warns COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 30 -t tracks/S-76D/197421AC.csv --wind 10kt,270 --stability II --gridstep 1640 -o integrated/wind_stability_class.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(wind_stability_class_warns PROPERTIES PASS_REGULAR_EXPRESSION "The stability classes are unvalidated")
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
# --terrain (src/terrain.hpp, terrain_ray in src/propagation.hpp).  terrain_files writes two small
# synthetic grids around the demo track's frame (origin 39.81043 -74.42631, +x along 210 degrees,
# its z = 0 at altitude 54.2651 ft): one level, at the height --terrain_offset -1250 and
# --ground_elevation 54.2651 take back to z = 0 exactly, and one of hills falling toward the
# approach's low end.  Level terrain must reproduce the run without it (point receivers, plates, and
# scattered observers stood on it), byte for byte in fact; over the hills the peak search must still
# equal the march, with the point receiver's ground model evaluated at each path's own height.
set(TERRAIN_DEMO -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --origin 39.81043 -74.42631 --frame_bearing 210 -L warn)
set(TERRAIN_PLATE --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc)
set(TERRAIN_FLAT --terrain integrated/terrain_flat.nc --terrain_offset -1250 --ground_elevation 54.2651)
set(TERRAIN_HILLS --terrain integrated/terrain_hills.nc --ground_elevation 54.2651)
add_test(NAME terrain_files COMMAND niceops_unit_tests "[.terrain_files]" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_files PROPERTIES FIXTURES_SETUP terrain_file)
add_test(NAME terrain_flat_point_reference COMMAND niceops ${TERRAIN_DEMO} --ground_effect --gridstep 492 -o integrated/terrain_flat_point_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_point COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_FLAT} --ground_effect --gridstep 492 -o integrated/terrain_flat_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_plate_reference COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_PLATE} --gridstep 492 -o integrated/terrain_flat_plate_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_plate COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_FLAT} ${TERRAIN_PLATE} --gridstep 492 -o integrated/terrain_flat_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_scattered_reference COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_PLATE} -b tests/data/coyle_approach_level.csv -o integrated/terrain_flat_scattered_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_scattered COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_FLAT} ${TERRAIN_PLATE} -b arrays/coyle_approach.csv --observers_on_terrain -o integrated/terrain_flat_scattered.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_flat_point_reference terrain_flat_plate_reference terrain_flat_scattered_reference PROPERTIES FIXTURES_SETUP terrain_flat_runs)
set_tests_properties(terrain_flat_point terrain_flat_plate terrain_flat_scattered PROPERTIES FIXTURES_REQUIRED terrain_file FIXTURES_SETUP terrain_flat_runs)
add_test(NAME terrain_flat_point_is_flat COMMAND golden_compare integrated/terrain_flat_point.csv integrated/terrain_flat_point_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_plate_is_flat COMMAND golden_compare integrated/terrain_flat_plate.csv integrated/terrain_flat_plate_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_scattered_is_flat COMMAND golden_compare integrated/terrain_flat_scattered.csv integrated/terrain_flat_scattered_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_flat_point_is_flat terrain_flat_plate_is_flat terrain_flat_scattered_is_flat PROPERTIES FIXTURES_REQUIRED terrain_flat_runs)
add_test(NAME terrain_hills_point COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} --ground_effect --gridstep 492 -o integrated/terrain_hills_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_point_marched COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} --ground_effect --gridstep 492 -D integrated/terrain_hills_point_dnl.csv -o integrated/terrain_hills_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_plate COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} ${TERRAIN_PLATE} --gridstep 492 -o integrated/terrain_hills_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_plate_marched COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} ${TERRAIN_PLATE} --gridstep 492 -D integrated/terrain_hills_plate_dnl.csv -o integrated/terrain_hills_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_hills_point terrain_hills_point_marched terrain_hills_plate terrain_hills_plate_marched PROPERTIES FIXTURES_REQUIRED terrain_file FIXTURES_SETUP terrain_hills_runs)
add_test(NAME terrain_hills_point_windowed_is_exact COMMAND golden_compare integrated/terrain_hills_point.csv integrated/terrain_hills_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_plate_windowed_is_exact COMMAND golden_compare integrated/terrain_hills_plate.csv integrated/terrain_hills_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_hills_point_windowed_is_exact terrain_hills_plate_windowed_is_exact PROPERTIES FIXTURES_REQUIRED terrain_hills_runs)
# ... and both are pinned, since the windowed/marched pair compares one binary with itself: passing
# the plate's receiver as a point receiver at its call site (the plate then measures the source from
# the line, not from its microphone) moved cells by 0.01 dB with every other terrain test green.
# Made by the build that added them (adopt/terrain); regenerate only with a deliberate change.
add_test(NAME terrain_hills_point_golden COMMAND golden_compare integrated/terrain_hills_point.csv golden/terrain_hills_point.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_plate_golden COMMAND golden_compare integrated/terrain_hills_plate.csv golden/terrain_hills_plate.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_hills_point_golden terrain_hills_plate_golden PROPERTIES FIXTURES_REQUIRED terrain_hills_runs)
# ... and the hills reach the ground model, not only the observers' heights: the same microphones
# stood on the hills through an observer file, with no --terrain, differ (golden_compare's mismatch
# is the pass) -- and without a ground model, where the terrain only places the observers, they agree.
add_test(NAME terrain_hills_scattered COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} --ground_effect -b arrays/coyle_approach.csv --observers_on_terrain -o integrated/terrain_hills_scattered.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_scattered_placed COMMAND niceops ${TERRAIN_DEMO} --ground_effect -b integrated/terrain_hills_observers.csv -o integrated/terrain_hills_scattered_placed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_bare COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} --no-ground_effect -b arrays/coyle_approach.csv --observers_on_terrain -o integrated/terrain_hills_bare.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_bare_placed COMMAND niceops ${TERRAIN_DEMO} --no-ground_effect -b integrated/terrain_hills_observers.csv -o integrated/terrain_hills_bare_placed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_hills_scattered terrain_hills_scattered_placed terrain_hills_bare terrain_hills_bare_placed PROPERTIES FIXTURES_REQUIRED terrain_file FIXTURES_SETUP terrain_placed_runs)
add_test(NAME terrain_reaches_the_ground_model COMMAND golden_compare integrated/terrain_hills_scattered.csv integrated/terrain_hills_scattered_placed.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_reaches_the_ground_model PROPERTIES FIXTURES_REQUIRED terrain_placed_runs PASS_REGULAR_EXPRESSION "mismatch at value")
add_test(NAME terrain_without_ground_places_observers COMMAND golden_compare integrated/terrain_hills_bare.csv integrated/terrain_hills_bare_placed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_without_ground_places_observers PROPERTIES FIXTURES_REQUIRED terrain_placed_runs)
# Refusals, passed on their message (a usage error or a missing file cannot satisfy them).
add_test(NAME terrain_needs_a_frame COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --terrain integrated/terrain_flat.nc -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_needs_a_frame PROPERTIES FIXTURES_REQUIRED terrain_file PASS_REGULAR_EXPRESSION "--terrain needs to know where")
add_test(NAME terrain_checks_ground_elevation COMMAND niceops ${TERRAIN_DEMO} --terrain integrated/terrain_flat.nc --terrain_offset -1250 --gridstep 984 -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_checks_ground_elevation PROPERTIES FIXTURES_REQUIRED terrain_file PASS_REGULAR_EXPRESSION "puts its z = 0 at altitude 54.2651")
add_test(NAME terrain_refuses_states_underground COMMAND niceops ${TERRAIN_DEMO} --terrain integrated/terrain_flat.nc --ground_elevation 54.2651 --gridstep 984 -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_refuses_states_underground PROPERTIES FIXTURES_REQUIRED terrain_file PASS_REGULAR_EXPRESSION "different vertical datums")
add_test(NAME terrain_checks_the_frame COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --origin 39.81043 -74.42631 --frame_bearing 200 ${TERRAIN_FLAT} --gridstep 984 -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_checks_the_frame PROPERTIES FIXTURES_REQUIRED terrain_file PASS_REGULAR_EXPRESSION "rms from its x/y")
add_test(NAME terrain_must_cover_the_observers COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --origin 39.95 -74.42631 --frame_bearing 210 ${TERRAIN_FLAT} -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_must_cover_the_observers PROPERTIES FIXTURES_REQUIRED terrain_file PASS_REGULAR_EXPRESSION "has no height at")
add_test(NAME terrain_refuses_refraction COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --origin 39.81043 -74.42631 --frame_bearing 210 ${TERRAIN_FLAT} --atmosphere atmospheres/B407_uniform.csv --refraction linear -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_refuses_refraction PROPERTIES FIXTURES_REQUIRED terrain_file PASS_REGULAR_EXPRESSION "does not combine with --refraction")
add_test(NAME terrain_offset_needs_terrain COMMAND niceops ${TERRAIN_DEMO} --terrain_offset -10 -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_offset_needs_terrain PROPERTIES PASS_REGULAR_EXPRESSION "need --terrain")
add_test(NAME terrain_warns_of_observer_disagreement COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_FLAT} ${TERRAIN_PLATE} -b arrays/coyle_approach.csv -o integrated/terrain_disagreement.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_warns_of_observer_disagreement PROPERTIES FIXTURES_REQUIRED terrain_file PASS_REGULAR_EXPRESSION "disagree by more than")
