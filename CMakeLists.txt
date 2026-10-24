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
# Turning it on is worth 18% of the loop and leaves the golden comparison
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

# GPU backend (src/gpu).  AUTO takes Metal on macOS with Clang and none
# elsewhere (GCC cannot build the Objective-C++ that drives Metal).  The
# kernels are compiled at run time from source embedded in the binary, so no
# shader toolchain is needed.  Where there is a device the tools that have a
# GPU path take it by default; --no-gpu keeps the CPU (docs/internals/gpu.md).
set(NICEOPS_GPU "AUTO" CACHE STRING "GPU backend: AUTO (Metal on macOS with Clang, else none), METAL or OFF")
set_property(CACHE NICEOPS_GPU PROPERTY STRINGS AUTO METAL OFF)
if(NICEOPS_GPU STREQUAL "AUTO")
  if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(NICEOPS_GPU_BACKEND METAL)
  else()
    set(NICEOPS_GPU_BACKEND OFF)
  endif()
else()
  set(NICEOPS_GPU_BACKEND ${NICEOPS_GPU})
endif()
if(NICEOPS_GPU_BACKEND STREQUAL "METAL")
  if(NOT APPLE)
    message(FATAL_ERROR "NICEOPS_GPU=METAL needs macOS")
  endif()
  # The backend is Objective-C++, which GCC does not build, and it shares
  # std::string and std::shared_ptr with the rest across gpu.hpp, so the
  # C++ compiler must be Clang too.
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    message(FATAL_ERROR "NICEOPS_GPU=METAL needs Clang as the C++ compiler, not ${CMAKE_CXX_COMPILER_ID}; "
                        "configure with -DNICEOPS_GPU=OFF to build without a GPU backend")
  endif()
  enable_language(OBJCXX)
  set(NICEOPS_GPU_SOURCE src/gpu/gpu_metal.mm)
  set_source_files_properties(src/gpu/gpu_metal.mm PROPERTIES COMPILE_OPTIONS "-fobjc-arc")
elseif(NICEOPS_GPU_BACKEND STREQUAL "OFF")
  set(NICEOPS_GPU_SOURCE src/gpu/gpu_disabled.cpp)
else()
  message(FATAL_ERROR "NICEOPS_GPU must be AUTO, METAL or OFF, not ${NICEOPS_GPU}")
endif()
message(STATUS "GPU backend: ${NICEOPS_GPU_BACKEND}")
# The kernel sources, embedded as strings (cmake/embed_kernels.cmake).  A new
# kernel file goes in this list.
set(NICEOPS_KERNEL_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/kernels/prelude_metal.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/kernels/common.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/kernels/band_energy_types.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/kernels/band_energy.kern
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/kernels/wedge.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/kernels/wedge_check.kern
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/kernels/point_rays_types.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/kernels/point_rays.kern
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/kernels/canyon_wave_types.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/kernels/canyon_wave.kern
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/kernels/fdtd_types.h
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/kernels/fdtd.kern
)

# Core library: noise database, trajectory, smoothing spline, plotting facade
add_library(niceops_core STATIC
  src/sphere_database.cpp
  src/atmosphere.cpp
  src/atmosphere_profile.cpp
  src/atmosphere_builder.cpp
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
  src/ground_classes.cpp
  src/plate_receiver.cpp
  src/geodetic.cpp
  src/terrain.cpp
  src/screening.cpp
  src/urban.cpp
  src/canyon_wave.cpp
  src/fdtd.cpp
  src/wedge.cpp
  src/double_wedge.cpp
  src/profile_pe.cpp
  src/edge_source.cpp
  src/buildings.cpp
  src/building_paths.cpp
  src/diffuse_facades.cpp
  src/geojson_buildings.cpp
  src/point_paths.cpp
  src/facade_relief.cpp
  src/grating_paths.cpp
  src/urban_point.cpp
  src/turbulence_scattering.cpp
  src/source_coherence.cpp
  src/edge_visibility.cpp
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
  src/contour.cpp
  src/raster.cpp
  src/sphere_energy.cpp
  src/state_condition.cpp
  src/procedure.cpp
  src/operations.cpp
  ${NICEOPS_GPU_SOURCE}
  src/gpu/band_energy.cpp
  src/gpu/point_rays.cpp
  src/gpu/band_stream.cpp
  ${CMAKE_CURRENT_BINARY_DIR}/generated/gpu_kernel_sources.cpp
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
# The list goes through as one argument, joined with '|': a ';' would split it.
string(REPLACE ";" "|" NICEOPS_KERNEL_SOURCES_JOINED "${NICEOPS_KERNEL_SOURCES}")
add_custom_command(
  OUTPUT ${NICEOPS_GENERATED_DIR}/gpu_kernel_sources.cpp
  COMMAND ${CMAKE_COMMAND} "-DSOURCES=${NICEOPS_KERNEL_SOURCES_JOINED}"
          -DOUTPUT=${NICEOPS_GENERATED_DIR}/gpu_kernel_sources.cpp
          -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_kernels.cmake
  DEPENDS ${NICEOPS_KERNEL_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/embed_kernels.cmake
  VERBATIM)
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
if(NICEOPS_GPU_BACKEND STREQUAL "METAL")
  target_link_libraries(niceops_core PUBLIC "-framework Metal" "-framework Foundation")
endif()
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
  tests/test_contour.cpp
  tests/test_procedure.cpp
  tests/test_iota.cpp
  tests/test_sphere_database.cpp
  tests/test_geodetic.cpp
  tests/test_terrain.cpp
  tests/test_screening.cpp
  tests/test_urban.cpp
  tests/test_canyon_wave.cpp
  tests/test_wedge.cpp
  tests/test_edge_source.cpp
  tests/test_fdtd.cpp
  tests/test_buildings.cpp
  tests/test_building_paths.cpp
  tests/test_diffuse_facades.cpp
  tests/test_geojson_buildings.cpp
  tests/test_point_paths.cpp
  tests/test_profile_pe.cpp
  tests/test_source_coherence.cpp
  tests/test_urban_point.cpp
  tests/test_turbulence_scattering.cpp
  tests/test_facade_relief.cpp
  tests/test_pitched_roofs.cpp
  tests/test_microphone_surfaces.cpp
  tests/test_gpu.cpp
  tests/test_edge_visibility.cpp
  tests/test_atmosphere.cpp
  tests/test_atmosphere_profile.cpp
  tests/test_atmosphere_builder.cpp
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
  tests/test_operations.cpp
  tests/test_ground_classes.cpp
  tests/test_legacy_repair.cpp
  tests/test_smoothing_spline.cpp
)
target_link_libraries(niceops_unit_tests niceops_core Catch2::Catch2WithMain)

# Golden output comparison tool (numeric, elementwise)
add_executable(golden_compare tests/golden_compare.cpp)
# What --epnl_file and --band_frames wrote, checked against each other and the SEL.
add_executable(check_band_frames tests/check_band_frames.cpp)
target_link_libraries(check_band_frames niceops_core)
# What --raster_file wrote, against the grid files the same run wrote.
add_executable(check_raster tests/check_raster.cpp)
target_link_libraries(check_raster niceops_core)
# What --provenance_spheres wrote, against the --provenance_file beside it.
add_executable(check_sphere_energy tests/check_sphere_energy.cpp)
target_link_libraries(check_sphere_energy niceops_core)
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
add_executable(niceops_fdtd_field tools/fdtd_field.cpp)
target_link_libraries(niceops_fdtd_field niceops_core CLI11::CLI11)
add_executable(niceops_procedure tools/procedure.cpp)
target_link_libraries(niceops_procedure niceops_core CLI11::CLI11)
add_executable(niceops_propagation_table tools/propagation_table.cpp)
target_link_libraries(niceops_propagation_table niceops_core CLI11::CLI11)
# An --atmosphere profile from measured columns, and a site roughness length
# from wind profiles (tools/build_atmosphere.cpp).
add_executable(niceops_build_atmosphere tools/build_atmosphere.cpp)
target_link_libraries(niceops_build_atmosphere niceops_core CLI11::CLI11)

add_executable(niceops_urban_table tools/urban_table.cpp)
target_link_libraries(niceops_urban_table niceops_core CLI11::CLI11)
add_executable(niceops_point_field tools/point_field.cpp)
target_link_libraries(niceops_point_field niceops_core CLI11::CLI11)
# The Attenborough et al. (1995) benchmark cases by the PE and the exact
# solutions (tools/benchmark.cpp); tools/plot_benchmark.py draws them.
add_executable(niceops_benchmark tools/benchmark.cpp)
target_link_libraries(niceops_benchmark niceops_core CLI11::CLI11 Threads::Threads)

# Add testing with CTest
enable_testing()
# CI builds with plotting only, so the stub the -DNICEOPS_PLOTTING=OFF build links is compiled here
# too (never linked: it would clash with plotting.cpp) and held to define every plt:: function
# plotting.cpp does.  Without it, a new entry point breaks only the minimal build's link.
if(NICEOPS_PLOTTING)
  add_library(niceops_plotting_stub OBJECT src/plotting_disabled.cpp)
  target_link_libraries(niceops_plotting_stub PRIVATE niceops_core)
  add_test(NAME plotting_stub_is_complete COMMAND ${CMAKE_COMMAND}
    -DNM=${CMAKE_NM}
    -DCORE=$<TARGET_FILE:niceops_core>
    -DSTUB=$<TARGET_OBJECTS:niceops_plotting_stub>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_plotting_stub.cmake)
endif()
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
# The demo case on a coarse grid: the reference for the runs below that show an option changes
# nothing (or something).  An identity holds or breaks on any grid, and the full one is
# regression_golden's; at 1640 ft a run costs a fraction of the default grid's.
set(DEMO_COARSE --gridstep 1640)
add_test(NAME demo_coarse COMMAND niceops --no-spectral --no-ground_effect ${DEMO_COARSE} -L warn -o integrated/demo_coarse.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(demo_coarse PROPERTIES FIXTURES_SETUP demo_coarse_run)
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
add_test(NAME screen_database_victim COMMAND ${CMAKE_COMMAND} -E copy databases/S-76D_fixed_load_factor.nod integrated/screen_victim.nod WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(screen_database_victim PROPERTIES FIXTURES_SETUP screen_victim)
add_test(NAME screen_database_refuses_overwriting_it COMMAND niceops -d integrated/screen_victim.nod --screen_database integrated/screen_victim.nod WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(screen_database_refuses_overwriting_it PROPERTIES FIXTURES_REQUIRED screen_victim PASS_REGULAR_EXPRESSION "would overwrite the database it reads")
add_test(NAME mask_suspect_cells_clean COMMAND niceops --no-spectral --no-ground_effect --mask_suspect_cells ${DEMO_COARSE} -L warn -o integrated/mask_suspect_clean.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(mask_suspect_cells_clean PROPERTIES FIXTURES_SETUP mask_clean_run)
# golden_compare's offset: a file against itself less 1 dB must mismatch, and the axes are exempt.
add_test(NAME golden_compare_offset_rejects COMMAND golden_compare integrated/197421_LAE.csv integrated/197421_LAE.csv 0 1 2 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(golden_compare_offset_rejects PROPERTIES FIXTURES_REQUIRED demo_run PASS_REGULAR_EXPRESSION "mismatch at value 240:")
add_test(NAME mask_suspect_cells_clean_is_golden COMMAND golden_compare integrated/mask_suspect_clean.csv integrated/demo_coarse.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(mask_suspect_cells_clean_is_golden PROPERTIES FIXTURES_REQUIRED "mask_clean_run;demo_coarse_run")
add_test(NAME mask_suspect_cells_refuses_self_test COMMAND niceops --self_test --mask_suspect_cells WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(mask_suspect_cells_refuses_self_test PROPERTIES PASS_REGULAR_EXPRESSION "does not combine with --mask_suspect_cells")
add_test(NAME screen_database_refuses_two_databases COMMAND niceops -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197421AC.csv -d databases/S-76D_M3.nod -d databases/S-76D_fixed_load_factor.nod -o integrated/screen_two.csv --screen_database integrated/screen_two_list.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(screen_database_refuses_two_databases PROPERTIES PASS_REGULAR_EXPRESSION "run it once per database")
# The Catch2 cases in shards, so that `ctest -j` runs them side by side: as one process they were the
# suite's longest test.  Each shard is its own process, so no case may lean on another's leftovers.
set(UNIT_TEST_SHARDS 12)
math(EXPR last_shard "${UNIT_TEST_SHARDS} - 1")
foreach(shard RANGE ${last_shard})
  add_test(NAME unit_test_${shard} COMMAND niceops_unit_tests "~[extended]" --shard-count ${UNIT_TEST_SHARDS} --shard-index ${shard} WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
endforeach()
# The cases tagged [extended] (accuracy against reference solutions, wide sweeps) run as one test,
# which CI leaves out with the other extended tests (tests/extended.cmake).
add_test(NAME unit_test_extended COMMAND niceops_unit_tests "[extended]" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# --manifest: the run's record.  The run is a plain test (exit status), the record is read back.
add_test(NAME manifest_run COMMAND niceops --no-spectral --no-ground_effect --gridstep 2000 -L warn -o integrated/manifest_demo.csv --manifest integrated/manifest_demo.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(manifest_run PROPERTIES FIXTURES_SETUP manifest_demo_run)
add_test(NAME manifest_records_run COMMAND ${CMAKE_COMMAND} -E cat integrated/manifest_demo.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(manifest_records_run PROPERTIES FIXTURES_REQUIRED manifest_demo_run PASS_REGULAR_EXPRESSION "\"role\": \"database\",[^}]*\"path\": \"databases/S-76D_M3\\.nod\",[^}]*\"sha256\": \"[0-9a-f]+\".*\"outputs\": \\[[^}]*\"path\": \"integrated/manifest_demo\\.csv\",[^}]*\"sha256\": \"[0-9a-f]+\".*\"exit_status\": 0")
# A run that stops still leaves its record, with the status it returned.
add_test(NAME manifest_records_failure COMMAND niceops -d databases/S-76D_fixed_load_factor.nod --self_test --mask_suspect_cells --manifest integrated/manifest_failed.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
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
add_test(NAME config_demo_is_golden COMMAND golden_compare integrated/config_demo.csv integrated/demo_coarse.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(config_demo_is_golden PROPERTIES FIXTURES_REQUIRED "config_demo_run;demo_coarse_run")
add_test(NAME config_refuses_an_unknown_key COMMAND niceops --config tests/data/misspelled_options.toml -o integrated/config_misspelled.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(config_refuses_an_unknown_key PROPERTIES PASS_REGULAR_EXPRESSION "gridstpe")
# --provenance_file: a load-factor ceiling of 1.05 g holds 600 of the demo's 6551 states on B407.
# It needs every state, so the run marches, and the SEL must be the windowed run's.
add_test(NAME provenance_run COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 --gridstep 2000 --max_load_factor 1.05 -L error -o integrated/provenance_sel.csv --provenance_file integrated/provenance.csv --provenance_spheres integrated/provenance_spheres.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME provenance_windowed COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 --gridstep 2000 --max_load_factor 1.05 -L error -o integrated/provenance_windowed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(provenance_run provenance_windowed PROPERTIES FIXTURES_SETUP provenance_runs)
add_test(NAME provenance_leaves_sel_alone COMMAND golden_compare integrated/provenance_sel.csv integrated/provenance_windowed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME provenance_shares_held_states COMMAND ${CMAKE_COMMAND} -E cat integrated/provenance.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(provenance_leaves_sel_alone provenance_shares_held_states PROPERTIES FIXTURES_REQUIRED provenance_runs)
set_tests_properties(provenance_shares_held_states PROPERTIES PASS_REGULAR_EXPRESSION "exposure_db,share_mu_clamped,.*\n-4921\\.[0-9]+,-4921\\.[0-9]+,62\\.583[0-9]*,0,0,0,0\\.0955[0-9]*,0,31,0\\.485")
# --provenance_spheres: every sphere's share, whose largest is the CSV's top sphere.
add_test(NAME provenance_spheres_match COMMAND check_sphere_energy integrated/provenance_spheres.nc integrated/provenance.csv 0.3048 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(provenance_spheres_match PROPERTIES FIXTURES_REQUIRED provenance_runs PASS_REGULAR_EXPRESSION "check_sphere_energy: [1-9][0-9]* of [0-9]+ observers match")
# Culled observers (--db_cutoff) heard nothing: a row of NaN shares and an exposure of -inf.
add_test(NAME provenance_spheres_silent_run COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 --gridstep 2000 --no-spectral --db_cutoff 80 -L error -o integrated/provenance_silent_sel.csv --provenance_file integrated/provenance_silent.csv --provenance_spheres integrated/provenance_silent.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(provenance_spheres_silent_run PROPERTIES FIXTURES_SETUP provenance_runs)
add_test(NAME provenance_spheres_silent COMMAND check_sphere_energy integrated/provenance_silent.nc integrated/provenance_silent.csv 0.3048 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(provenance_spheres_silent PROPERTIES FIXTURES_REQUIRED provenance_runs PASS_REGULAR_EXPRESSION "check_sphere_energy: [1-9][0-9]* of 25 observers match, [1-9][0-9]* heard nothing")
add_test(NAME provenance_spheres_need_the_file COMMAND niceops --provenance_spheres integrated/provenance_alone.nc -o integrated/provenance_alone_sel.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(provenance_spheres_need_the_file PROPERTIES PASS_REGULAR_EXPRESSION "--provenance_spheres writes beside --provenance_file: give that too")
add_test(NAME provenance_refuses_ambient COMMAND niceops -a 40 --provenance_file integrated/provenance_ambient.csv -o integrated/provenance_ambient_sel.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(provenance_refuses_ambient PROPERTIES PASS_REGULAR_EXPRESSION "does not combine with -a")
add_test(NAME provenance_refuses_two_databases COMMAND niceops -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197421AC.csv -d databases/S-76D_M3.nod -d databases/S-76D_fixed_load_factor.nod -o integrated/provenance_two.csv --provenance_file integrated/provenance_two_list.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(provenance_refuses_two_databases PROPERTIES PASS_REGULAR_EXPRESSION "run it once per database")
# The database cache: each distinct database loads once, however its trajectories are interleaved
# (A B A B: the front one is held throughout, the other until its last trajectory) and however its
# path is spelled.  Every load opens the file ("Opening netCDF file"), and info lines
# reach only the log, so each run writes a fresh log of its own and the check reads it back.
add_test(NAME database_cache_run COMMAND sh -c "rm -f integrated/database_cache.log && $<TARGET_FILE:niceops> --no-spectral --no-ground_effect --gridstep 2000 -L info -l integrated/database_cache.log -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197422AC.csv -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197422AC.csv -d databases/S-76D_M3.nod -d databases/B407_fixed_load_factor.nod -d databases/S-76D_M3.nod -d databases/B407_fixed_load_factor.nod -v databases/S-76D.json -v databases/B407.json -v databases/S-76D.json -v databases/B407.json -o integrated/database_cache_a1.csv -o integrated/database_cache_b1.csv -o integrated/database_cache_a2.csv -o integrated/database_cache_b2.csv" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(database_cache_run PROPERTIES FIXTURES_SETUP database_cache_run)
add_test(NAME database_cache_loads_each_once COMMAND ${CMAKE_COMMAND} -E cat integrated/database_cache.log WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(database_cache_loads_each_once PROPERTIES FIXTURES_REQUIRED database_cache_run PASS_REGULAR_EXPRESSION "Loading database databases/B407_fixed_load_factor\\.nod for tracks/S-76D/197422AC\\.csv" FAIL_REGULAR_EXPRESSION "Opening netCDF file: databases/S-76D_M3\\.nod.*Opening netCDF file: databases/S-76D_M3\\.nod;Opening netCDF file: databases/B407_fixed_load_factor\\.nod.*Opening netCDF file: databases/B407_fixed_load_factor\\.nod")
# A held database gives a trajectory what it gave the same trajectory before.
add_test(NAME database_cache_reuse_is_identical COMMAND ${CMAKE_COMMAND} -E compare_files integrated/database_cache_a1.csv integrated/database_cache_a2.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(database_cache_reuse_is_identical PROPERTIES FIXTURES_REQUIRED database_cache_run)
add_test(NAME database_cache_reuse_is_identical_second COMMAND ${CMAKE_COMMAND} -E compare_files integrated/database_cache_b1.csv integrated/database_cache_b2.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(database_cache_reuse_is_identical_second PROPERTIES FIXTURES_REQUIRED database_cache_run)
# Two spellings of one database are one database: loaded once, and the provenance takes them.
add_test(NAME database_spellings_run COMMAND sh -c "rm -f integrated/database_spellings.log && $<TARGET_FILE:niceops> --no-spectral --no-ground_effect --gridstep 2000 -L info -l integrated/database_spellings.log -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197422AC.csv -d databases/S-76D_M3.nod -d databases/../databases/S-76D_M3.nod -o integrated/database_spellings.csv --provenance_file integrated/database_spellings_provenance.csv" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(database_spellings_run PROPERTIES FIXTURES_SETUP database_spellings_run)
add_test(NAME database_spellings_load_once COMMAND ${CMAKE_COMMAND} -E cat integrated/database_spellings.log WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(database_spellings_load_once PROPERTIES FIXTURES_REQUIRED database_spellings_run PASS_REGULAR_EXPRESSION "Opening netCDF file: databases/S-76D_M3\\.nod" FAIL_REGULAR_EXPRESSION "Loading database;Opening netCDF file: databases/\\.\\./")
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
add_test(NAME ambient_negligible_is_the_golden COMMAND niceops --no-spectral --no-ground_effect -a -200 ${DEMO_COARSE} -o integrated/ambient_negligible.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ambient_negligible_matches COMMAND golden_compare integrated/ambient_negligible.csv integrated/demo_coarse.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ambient_negligible_is_the_golden PROPERTIES FIXTURES_SETUP ambient_run)
set_tests_properties(ambient_negligible_matches PROPERTIES FIXTURES_REQUIRED "ambient_run;demo_coarse_run")
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
add_test(NAME period_night COMMAND niceops ${PERIOD_RUN} --start_utc 3600 --lden_file integrated/period_night_lden.csv --lnight_file integrated/period_night_lnight.csv --laeq_file integrated/period_night_laeq.csv --laeq_window 23 7 --lday_file integrated/period_night_lday.csv -o integrated/period_night.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME period_evening COMMAND niceops ${PERIOD_RUN} --start_utc 72000 --lden_file integrated/period_evening_lden.csv --levening_file integrated/period_evening_levening.csv -o integrated/period_evening.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(period_noon period_noon_day_window period_night period_evening PROPERTIES FIXTURES_SETUP period_runs)
# Each over its own length, against Lden over the day: a flight in one period gives
# Lday = Lden + 10 lg 2, Levening = Lden - 5 + 10 lg 6, Lnight = Lden - 10 + 10 lg 3.
add_test(NAME period_lday_over_12_hours COMMAND golden_compare integrated/period_noon_lday.csv integrated/period_noon_lden.csv 0 3.010299956639812 2 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME period_levening_over_4_hours COMMAND golden_compare integrated/period_evening_levening.csv integrated/period_evening_lden.csv 0 2.781512503836437 2 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME period_lnight_over_8_hours COMMAND golden_compare integrated/period_night_lnight.csv integrated/period_night_lden.csv 0 -5.228787452803376 2 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(period_lday_over_12_hours period_levening_over_4_hours period_lnight_over_8_hours PROPERTIES FIXTURES_REQUIRED period_runs)
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
# --contour_file: GeoJSON polygons of the footprint at the --cmin..--cmax levels, on the geodetic
# demo track (whose frame the run knows).  The polygons' geometry is pinned by the unit tests
# ([contour]); here the file, its levels and the refusals.
add_test(NAME contour_run COMMAND niceops ${PERIOD_RUN} -o integrated/contour_sel.csv --contour_file integrated/contours.geojson WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(contour_run PROPERTIES FIXTURES_SETUP contour_runs)
add_test(NAME contour_writes_geojson COMMAND ${CMAKE_COMMAND} -E cat integrated/contours.geojson WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(contour_writes_geojson PROPERTIES FIXTURES_REQUIRED contour_runs PASS_REGULAR_EXPRESSION "^\\{\"type\":\"FeatureCollection\",\"features\":\\[\\{\"type\":\"Feature\",\"properties\":\\{\"metric\":\"sel\",\"level_db\":65\\.0,\"area_m2\":[0-9.e+]+,.*\"geometry\":\\{\"type\":\"MultiPolygon\",\"coordinates\":\\[\\[\\[\\[-74\\.[0-9]+,39\\.[0-9]+\\]")
add_test(NAME contour_dnl COMMAND niceops ${PERIOD_RUN} --start_utc 43200 -D integrated/contour_dnl_levels.csv --contour_metric dnl --cmin 20 --cmax 30 --cstep 5 -o integrated/contour_dnl_sel.csv --contour_file integrated/contours_dnl.geojson WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(contour_dnl PROPERTIES FIXTURES_SETUP contour_runs)
add_test(NAME contour_dnl_levels COMMAND ${CMAKE_COMMAND} -E cat integrated/contours_dnl.geojson WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(contour_dnl_levels PROPERTIES FIXTURES_REQUIRED contour_runs PASS_REGULAR_EXPRESSION "\"metric\":\"dnl\",\"level_db\":20\\.0,\"area_m2\":[1-9].*\"level_db\":25\\.0.*\"level_db\":30\\.0")
add_test(NAME contour_needs_a_frame COMMAND niceops --gridstep 1640 -o integrated/contour_noframe.csv --contour_file integrated/contours_noframe.geojson WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(contour_needs_a_frame PROPERTIES PASS_REGULAR_EXPRESSION "--contour_file writes WGS84 polygons and needs the run's frame")
add_test(NAME contour_needs_the_grid COMMAND niceops ${PERIOD_RUN} -b tests/data/mixed_flush.csv -o integrated/contour_scattered.csv --contour_file integrated/contours_scattered.geojson WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(contour_needs_the_grid PROPERTIES PASS_REGULAR_EXPRESSION "scattered observers \\(-b\\) have none")
add_test(NAME contour_metric_needs_its_output COMMAND niceops ${PERIOD_RUN} --contour_metric lamax -o integrated/contour_nolamax.csv --contour_file integrated/contours_nolamax.geojson WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(contour_metric_needs_its_output PROPERTIES PASS_REGULAR_EXPRESSION "--contour_metric lamax contours that metric's output: ask for it too")
# --raster_file: every gridded level of the run in one CF netCDF, each field the grid file it
# mirrors (a night level at noon is all silence, written as the fill value), and the latitude and
# longitude the --latlon companion's.  -m's last column is the area at or above 65 dB.
add_test(NAME raster_clean COMMAND ${CMAKE_COMMAND} -E rm -f integrated/raster_metrics.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(raster_clean PROPERTIES FIXTURES_SETUP raster_clean)
add_test(NAME raster_run COMMAND niceops ${PERIOD_RUN} --start_utc 43200 --lday_file integrated/raster_lday.csv --lnight_file integrated/raster_lnight.csv --lamax_file integrated/raster_lamax.csv --lamax_weighting slow --ta_file integrated/raster_ta.csv --ta_threshold 65 --grid -4921.26 4921.26 -3281.26 3281.26 --latlon -o integrated/raster_sel.csv --raster_file integrated/raster.nc -m integrated/raster_metrics.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(raster_run PROPERTIES FIXTURES_SETUP raster_runs FIXTURES_REQUIRED raster_clean)
add_test(NAME raster_matches_grids COMMAND check_raster integrated/raster.nc 0.3048 --latlon integrated/raster_sel_latlon.csv sel=integrated/raster_sel.csv lday=integrated/raster_lday.csv lnight=integrated/raster_lnight.csv lamax=integrated/raster_lamax.csv time_above=integrated/raster_ta.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(raster_matches_grids PROPERTIES FIXTURES_REQUIRED raster_runs PASS_REGULAR_EXPRESSION "check_raster: 5 fields match")
add_test(NAME metrics_file_area COMMAND ${CMAKE_COMMAND} -E cat integrated/raster_metrics.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(metrics_file_area PROPERTIES FIXTURES_REQUIRED raster_runs PASS_REGULAR_EXPRESSION "^demo_geodetic_noclock,73\\.278[0-9]+,82\\.0017[0-9]+,0\\.62857[0-9]+,4073883[0-9]\\.[0-9]+\n?$")
add_test(NAME raster_without_frame COMMAND niceops --gridstep 1640 -o integrated/raster_noframe_sel.csv --raster_file integrated/raster_noframe.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(raster_without_frame PROPERTIES FIXTURES_SETUP raster_runs)
add_test(NAME raster_without_frame_is_local COMMAND check_raster integrated/raster_noframe.nc 0.3048 --no-latlon sel=integrated/raster_noframe_sel.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(raster_without_frame_is_local PROPERTIES FIXTURES_REQUIRED raster_runs PASS_REGULAR_EXPRESSION "check_raster: 1 fields match")
# A skipped trajectory writes no -o, so the raster holds no SEL (its vector is only zeros); the
# coverage share, written as its text file is, still goes in.
add_test(NAME raster_skipped_trajectory COMMAND niceops --gridstep 1640 -o integrated/raster_skipped_sel.csv --coverage_file integrated/raster_skipped_coverage.csv --raster_file integrated/raster_skipped.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(raster_skipped_trajectory PROPERTIES FIXTURES_SETUP raster_runs WILL_FAIL TRUE)
add_test(NAME raster_skipped_has_no_sel COMMAND check_raster integrated/raster_skipped.nc 0.3048 --no-latlon --absent sel unmeasured_share=integrated/raster_skipped_coverage.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(raster_skipped_has_no_sel PROPERTIES FIXTURES_REQUIRED raster_runs PASS_REGULAR_EXPRESSION "check_raster: 1 fields match")
add_test(NAME raster_needs_the_grid COMMAND niceops ${PERIOD_RUN} -b tests/data/mixed_flush.csv -o integrated/raster_scattered.csv --raster_file integrated/raster_scattered.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(raster_needs_the_grid PROPERTIES PASS_REGULAR_EXPRESSION "--raster_file writes the observer grid; scattered observers \\(-b\\) have none")
# The air per trajectory: two copies of the demo track, each with its own humidity (and its own
# --atmosphere), must each be the run of that track alone with that air; and the two differ.
set(AIR_RUN -d databases/B407_fixed_load_factor.nod -v databases/B407.json --gridstep 1640 -L warn)
add_test(NAME air_per_track COMMAND niceops ${AIR_RUN} -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197421AC.csv --humidity 20 --humidity 60 -o integrated/air_track_a.csv -o integrated/air_track_b.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_alone_20 COMMAND niceops ${AIR_RUN} -t tracks/S-76D/197421AC.csv --humidity 20 -o integrated/air_alone_20.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_alone_60 COMMAND niceops ${AIR_RUN} -t tracks/S-76D/197421AC.csv --humidity 60 -o integrated/air_alone_60.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_profile_per_track COMMAND niceops ${AIR_RUN} -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/isa.csv --atmosphere atmospheres/B407_uniform.csv -o integrated/air_profile_a.csv -o integrated/air_profile_b.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_profile_alone_isa COMMAND niceops ${AIR_RUN} -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/isa.csv -o integrated/air_profile_alone_isa.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_profile_alone_uniform COMMAND niceops ${AIR_RUN} -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/B407_uniform.csv -o integrated/air_profile_alone_uniform.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(air_per_track air_alone_20 air_alone_60 air_profile_per_track air_profile_alone_isa air_profile_alone_uniform PROPERTIES FIXTURES_SETUP air_runs)
add_test(NAME air_first_track_is_its_own COMMAND golden_compare integrated/air_track_a.csv integrated/air_alone_20.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_second_track_is_its_own COMMAND golden_compare integrated/air_track_b.csv integrated/air_alone_60.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_tracks_differ COMMAND golden_compare integrated/air_track_a.csv integrated/air_track_b.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_first_profile_is_its_own COMMAND golden_compare integrated/air_profile_a.csv integrated/air_profile_alone_isa.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_second_profile_is_its_own COMMAND golden_compare integrated/air_profile_b.csv integrated/air_profile_alone_uniform.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(air_first_track_is_its_own air_second_track_is_its_own air_tracks_differ air_first_profile_is_its_own air_second_profile_is_its_own PROPERTIES FIXTURES_REQUIRED air_runs)
set_tests_properties(air_tracks_differ PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value")
# Wind, temperature and pressure per track likewise (through --refraction, so the wind matters),
# and the clock: a noon copy and a 01:00 copy give the 01:00 copy's Lnight exactly, the noon one
# adding nothing.
set(AIR_WIND -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 50 --refraction linear --frame_bearing 210 --gridstep 1640 -L warn)
add_test(NAME air_wind_per_track COMMAND niceops ${AIR_WIND} -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197421AC.csv --wind 15kt,30 --wind 5kt,200 --temperature 280 --temperature 300 --pressure 95 --pressure 101 -o integrated/air_wind_a.csv -o integrated/air_wind_b.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_wind_alone_a COMMAND niceops ${AIR_WIND} -t tracks/S-76D/197421AC.csv --wind 15kt,30 --temperature 280 --pressure 95 -o integrated/air_wind_alone_a.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_wind_alone_b COMMAND niceops ${AIR_WIND} -t tracks/S-76D/197421AC.csv --wind 5kt,200 --temperature 300 --pressure 101 -o integrated/air_wind_alone_b.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_clock_per_track COMMAND niceops ${PERIOD_RUN} -t tests/data/demo_geodetic_noclock.csv --start_utc 43200 --start_utc 3600 --lnight_file integrated/air_clock_both.csv -o integrated/air_clock_a.csv -o integrated/air_clock_b.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_clock_alone COMMAND niceops ${PERIOD_RUN} --start_utc 3600 --lnight_file integrated/air_clock_alone.csv -o integrated/air_clock_alone_sel.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(air_wind_per_track air_wind_alone_a air_wind_alone_b air_clock_per_track air_clock_alone PROPERTIES FIXTURES_SETUP air_more_runs)
add_test(NAME air_first_wind_is_its_own COMMAND golden_compare integrated/air_wind_a.csv integrated/air_wind_alone_a.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_second_wind_is_its_own COMMAND golden_compare integrated/air_wind_b.csv integrated/air_wind_alone_b.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME air_clock_is_each_tracks_own COMMAND golden_compare integrated/air_clock_both.csv integrated/air_clock_alone.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(air_first_wind_is_its_own air_second_wind_is_its_own air_clock_is_each_tracks_own PROPERTIES FIXTURES_REQUIRED air_more_runs)
# The wind's bearing warning reaches a blowing wind that follows a calm one.
add_test(NAME air_wind_bearing_warns_after_calm COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 50 --refraction linear --gridstep 4000 -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197421AC.csv --wind 0kt,0 --wind 15kt,30 -o integrated/air_bearing_a.csv -o integrated/air_bearing_b.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(air_wind_bearing_warns_after_calm PROPERTIES PASS_REGULAR_EXPRESSION "the direction is a compass bearing")
add_test(NAME air_count_must_match COMMAND niceops ${AIR_RUN} -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197421AC.csv --humidity 20 --humidity 40 --humidity 60 -o integrated/air_count_a.csv -o integrated/air_count_b.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(air_count_must_match PROPERTIES PASS_REGULAR_EXPRESSION "Got 3 --humidity values for 2 trajectories")
# niceops_procedure: a descending, decelerating approach to a hover, flown into a WGS84 track
# with velocity and heading columns, which niceops then runs as it is.
add_test(NAME procedure_build COMMAND niceops_procedure tests/data/procedure_approach.json -o integrated/procedure_approach.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(procedure_build PROPERTIES FIXTURES_SETUP procedure_track)
add_test(NAME procedure_runs COMMAND niceops -t integrated/procedure_approach.csv --origin 39.8291 -74.4100 --frame_bearing 210 --gridstep 1640 -L warn -o integrated/procedure_approach_sel.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(procedure_runs PROPERTIES FIXTURES_REQUIRED procedure_track)
# niceops_procedure --check: 150 kt is past B407's advance ratios and a 50-degree bank past
# 1.3 g; the states it reports clamped and held, and when, are the ones niceops then counts.
add_test(NAME procedure_check COMMAND niceops_procedure tests/data/procedure_envelope.json -o integrated/procedure_envelope.csv --check databases/B407_fixed_load_factor.nod --vehicle databases/B407.json --max_load_factor 1.3 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(procedure_check PROPERTIES FIXTURES_SETUP procedure_envelope PASS_REGULAR_EXPRESSION "122 of 301 states \\(every 0\\.20 s\\) clamped or held\\.\n  advance ratio clamped to the database's range: 87 states, t = 0\\.0-17\\.2 s\n.*load factor held at --max_load_factor: 35 states, t = 44\\.0-50\\.8 s")
add_test(NAME procedure_check_agrees_with_niceops COMMAND niceops -t integrated/procedure_envelope.csv -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 --max_load_factor 1.3 --gridstep 4000 -L warn -o integrated/procedure_envelope_sel.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(procedure_check_agrees_with_niceops PROPERTIES FIXTURES_REQUIRED procedure_envelope PASS_REGULAR_EXPRESSION "advance ratio 87, flight-path angle 0, thrust coefficient 0 of 301 states.*35 of 301 states exceed 1\\.3 g")
add_test(NAME procedure_check_refuses_wrong_vehicle COMMAND niceops_procedure tests/data/procedure_envelope.json -o integrated/procedure_envelope_s76.csv --check databases/B407_fixed_load_factor.nod --vehicle databases/S-76D.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(procedure_check_refuses_wrong_vehicle PROPERTIES PASS_REGULAR_EXPRESSION "Vehicle main rotor tip speed disagrees: databases/S-76D.json gives 220 m/s.*so the run is refused")
add_test(NAME procedure_check_needs_vehicle COMMAND niceops_procedure tests/data/procedure_envelope.json -o integrated/procedure_envelope_nov.csv --check databases/B407_fixed_load_factor.nod WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(procedure_check_needs_vehicle PROPERTIES PASS_REGULAR_EXPRESSION "--check maps the states through the vehicle's rotor: give --vehicle")
add_test(NAME procedure_check_needs_a_positive_load_factor COMMAND niceops_procedure tests/data/procedure_envelope.json -o integrated/procedure_envelope_zero.csv --check databases/B407_fixed_load_factor.nod --vehicle databases/B407.json --max_load_factor 0 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(procedure_check_needs_a_positive_load_factor PROPERTIES PASS_REGULAR_EXPRESSION "--max_load_factor: Value 0 not in range")
add_test(NAME procedure_refuses_going_underground COMMAND niceops_procedure tests/data/procedure_underground.json -o integrated/procedure_underground.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(procedure_refuses_going_underground PROPERTIES PASS_REGULAR_EXPRESSION "segment 2 goes below the ground")
# A track below the observers' plane (here an altitude datum 5000 ft too high) ran silently.
add_test(NAME track_below_ground_warns COMMAND niceops -t tests/data/demo_geodetic_noclock.csv --origin 39.82914323 -74.40996141 --ground_elevation 5000 --gridstep 1640 -L warn -o integrated/below_ground.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(track_below_ground_warns PROPERTIES PASS_REGULAR_EXPRESSION "states of tests/data/demo_geodetic_noclock\\.csv are below the ground")
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
# A static hold has no clock of its own.  Its states read UTC midnight, so -D gave every hold
# the night penalty (+10 dB) without a word; it is refused now, and --start_utc dates it.
add_test(NAME static_hold_dnl_needs_a_clock COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tests/data/static_hover.csv --gridstep 1640 -L warn -D integrated/static_hover_noclock_dnl.csv -o integrated/static_hover_noclock.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(static_hold_dnl_needs_a_clock PROPERTIES PASS_REGULAR_EXPRESSION "-D needs the time of day")
add_test(NAME static_hold_dnl_with_start_utc COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tests/data/static_hover.csv --gridstep 1640 -L warn --start_utc 43200 -D integrated/static_hover_dated_dnl.csv -o integrated/static_hover_dated.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
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
# --operations.  The fixtures in tests/data/operations/ name their files relative to that folder.
# A count of 1 is the -t run bit for bit (on the coarse grid; the full one is regression_golden's).
add_test(NAME operations_count_one_run COMMAND niceops --no-spectral --no-ground_effect -d databases/S-76D_M3.nod -v databases/S-76D.json --operations tests/data/operations/count1.csv ${DEMO_COARSE} -o integrated/ops_count1.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_count_one_run PROPERTIES FIXTURES_SETUP ops_count1_run)
add_test(NAME operations_count_one_is_bit_identical COMMAND ${CMAKE_COMMAND} -E compare_files integrated/ops_count1.csv integrated/demo_coarse.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_count_one_is_bit_identical PROPERTIES FIXTURES_REQUIRED "ops_count1_run;demo_coarse_run")
# The run totals of one track flown once by -t, and by operations files of count 1, 10 and 3, and
# by three repeated -t.  Per observer a total starts at zero, so a count of 3 gives fl(3x), which
# is what 0 + x + x + x gives: every total but the provenance is bit-identical to repeated runs.
set(OPS_RUN -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 --gridstep 1640 -L warn)
foreach(case t c1 c10 c3 rep)
  set(OPS_TOTALS_${case} -D integrated/ops_${case}_dnl.csv --lden_file integrated/ops_${case}_lden.csv --cnel_file integrated/ops_${case}_cnel.csv --lamax_file integrated/ops_${case}_lamax.csv --nax_file integrated/ops_${case}_nax.csv --nax_threshold 45 --ta_file integrated/ops_${case}_ta.csv --ta_threshold 40)
endforeach()
add_test(NAME operations_totals_t COMMAND niceops ${OPS_RUN} -t tracks/S-76D/197421AC.csv ${OPS_TOTALS_t} -o integrated/ops_t.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_totals_count1 COMMAND niceops ${OPS_RUN} --operations tests/data/operations/count1.csv ${OPS_TOTALS_c1} -o integrated/ops_c1.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_totals_count10 COMMAND niceops ${OPS_RUN} --operations tests/data/operations/count10.csv ${OPS_TOTALS_c10} -o integrated/ops_c10.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_totals_count3 COMMAND niceops ${OPS_RUN} --operations tests/data/operations/count3.csv ${OPS_TOTALS_c3} -o integrated/ops_c3.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_totals_repeated_t COMMAND niceops ${OPS_RUN} -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197421AC.csv ${OPS_TOTALS_rep} -o integrated/ops_rep_a.csv -o integrated/ops_rep_b.csv -o integrated/ops_rep_c.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# The coverage and provenance pair, on a database with gated cells.
set(OPS_GATED -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 --gridstep 2000 -L error)
add_test(NAME operations_gated_count3 COMMAND niceops ${OPS_GATED} --operations tests/data/operations/count3.csv --coverage_file integrated/ops_c3_coverage.csv --provenance_file integrated/ops_c3_provenance.csv -o integrated/ops_c3_gated.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_gated_repeated_t COMMAND niceops ${OPS_GATED} -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197421AC.csv --coverage_file integrated/ops_rep_coverage.csv --provenance_file integrated/ops_rep_provenance.csv -o integrated/ops_rep_ga.csv -o integrated/ops_rep_gb.csv -o integrated/ops_rep_gc.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_totals_t operations_totals_count1 operations_totals_count10 operations_totals_count3 operations_totals_repeated_t operations_gated_count3 operations_gated_repeated_t PROPERTIES FIXTURES_SETUP ops_totals_runs)
add_test(NAME operations_count_one_totals_are_t COMMAND ${CMAKE_COMMAND} "-DPAIRS=integrated/ops_c1.csv;integrated/ops_t.csv;integrated/ops_c1_dnl.csv;integrated/ops_t_dnl.csv;integrated/ops_c1_lden.csv;integrated/ops_t_lden.csv;integrated/ops_c1_cnel.csv;integrated/ops_t_cnel.csv;integrated/ops_c1_nax.csv;integrated/ops_t_nax.csv;integrated/ops_c1_ta.csv;integrated/ops_t_ta.csv;integrated/ops_c1_lamax.csv;integrated/ops_t_lamax.csv" -P tests/compare_pairs.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# A count of 10 is +10 dB on the day metrics (at 1e-9: the multiply and the logarithm each
# round), and leaves the SEL file (one operation's) and LAmax as they were.
add_test(NAME operations_count_ten_adds_10_db_dnl COMMAND golden_compare integrated/ops_c10_dnl.csv integrated/ops_t_dnl.csv 0 10 2 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_count_ten_adds_10_db_lden COMMAND golden_compare integrated/ops_c10_lden.csv integrated/ops_t_lden.csv 0 10 2 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_count_ten_adds_10_db_cnel COMMAND golden_compare integrated/ops_c10_cnel.csv integrated/ops_t_cnel.csv 0 10 2 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_count_ten_leaves_sel_and_lamax COMMAND ${CMAKE_COMMAND} "-DPAIRS=integrated/ops_c10.csv;integrated/ops_t.csv;integrated/ops_c10_lamax.csv;integrated/ops_t_lamax.csv" -P tests/compare_pairs.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_count_three_equals_repeated_t COMMAND ${CMAKE_COMMAND} "-DPAIRS=integrated/ops_c3_dnl.csv;integrated/ops_rep_dnl.csv;integrated/ops_c3_lden.csv;integrated/ops_rep_lden.csv;integrated/ops_c3_cnel.csv;integrated/ops_rep_cnel.csv;integrated/ops_c3_nax.csv;integrated/ops_rep_nax.csv;integrated/ops_c3_ta.csv;integrated/ops_rep_ta.csv;integrated/ops_c3_lamax.csv;integrated/ops_rep_lamax.csv;integrated/ops_c3_coverage.csv;integrated/ops_rep_coverage.csv" -P tests/compare_pairs.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_count_three_provenance COMMAND golden_compare integrated/ops_c3_provenance.csv integrated/ops_rep_provenance.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_count_one_totals_are_t operations_count_ten_adds_10_db_dnl operations_count_ten_adds_10_db_lden operations_count_ten_adds_10_db_cnel operations_count_ten_leaves_sel_and_lamax operations_count_three_equals_repeated_t operations_count_three_provenance PROPERTIES FIXTURES_REQUIRED ops_totals_runs)
# local_time dates each operation in place of its track's clock: 23:30 is DNL's night (+10 dB
# on 12:00), 12:00 is the track's own daytime clock's DNL, and a static hold at 12:00 is the
# hold dated by --start_utc 43200 (UTC is local here).
add_test(NAME operations_local_noon COMMAND niceops ${OPS_RUN} --operations tests/data/operations/local_noon.csv -D integrated/ops_noon_dnl.csv -o integrated/ops_noon.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_local_night COMMAND niceops ${OPS_RUN} --operations tests/data/operations/local_night.csv -D integrated/ops_night_dnl.csv -o integrated/ops_night.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_local_hold COMMAND niceops ${OPS_RUN} --operations tests/data/operations/static_local.csv -D integrated/ops_hold_dnl.csv -o integrated/ops_hold.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_local_noon operations_local_night operations_local_hold PROPERTIES FIXTURES_SETUP ops_local_runs)
add_test(NAME operations_local_time_night_adds_10 COMMAND golden_compare integrated/ops_night_dnl.csv integrated/ops_noon_dnl.csv 0 10 2 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_local_time_night_adds_10 PROPERTIES FIXTURES_REQUIRED ops_local_runs)
add_test(NAME operations_local_time_replaces_clock COMMAND ${CMAKE_COMMAND} -E compare_files integrated/ops_noon_dnl.csv integrated/ops_t_dnl.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_local_time_replaces_clock PROPERTIES FIXTURES_REQUIRED "ops_local_runs;ops_totals_runs")
add_test(NAME operations_local_time_dates_a_static_hold COMMAND ${CMAKE_COMMAND} -E compare_files integrated/ops_hold_dnl.csv integrated/static_hover_dated_dnl.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(static_hold_dnl_with_start_utc PROPERTIES FIXTURES_SETUP static_hold_dated_run)
set_tests_properties(operations_local_time_dates_a_static_hold PROPERTIES FIXTURES_REQUIRED "ops_local_runs;static_hold_dated_run")
# Day, evening and night counts weight the plain SEL in the clock's place: a day count of 1 is the
# demo track flown by day byte for byte, a night count of 1 is it flown at 23:30 at 1e-9, and
# 2, 1, 1 is 1, 0, 0 plus 10 lg(2 + 10^0.5 + 10) on Lden, 10 lg 13 on DNL (evening counted as day)
# and 10 lg 15 on CNEL.  A track with no clock runs on counts.
set(OPS_PERIODS_night -D integrated/ops_pnight_dnl.csv --lden_file integrated/ops_pnight_lden.csv --cnel_file integrated/ops_pnight_cnel.csv)
set(OPS_PERIODS_mix -D integrated/ops_pmix_dnl.csv --lden_file integrated/ops_pmix_lden.csv --cnel_file integrated/ops_pmix_cnel.csv)
add_test(NAME operations_period_day_run COMMAND sh -c "rm -f integrated/ops_pday.log && $<TARGET_FILE:niceops> -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 --gridstep 1640 -L info -l integrated/ops_pday.log --operations tests/data/operations/period_day.csv -D integrated/ops_pday_dnl.csv --lden_file integrated/ops_pday_lden.csv --cnel_file integrated/ops_pday_cnel.csv -o integrated/ops_pday.csv" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_period_ignores_the_track_clock COMMAND ${CMAKE_COMMAND} -E cat integrated/ops_pday.log WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_period_ignores_the_track_clock PROPERTIES FIXTURES_REQUIRED ops_period_runs PASS_REGULAR_EXPRESSION "\\[info\\] 197421AC: its day, evening and night counts weight the day metrics, so the clock of tracks/S-76D/197421AC\\.csv is not read for them")
add_test(NAME operations_period_night_run COMMAND niceops ${OPS_RUN} --operations tests/data/operations/period_night.csv ${OPS_PERIODS_night} -o integrated/ops_pnight.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_period_mix_run COMMAND niceops ${OPS_RUN} --operations tests/data/operations/period_mix.csv ${OPS_PERIODS_mix} -o integrated/ops_pmix.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_period_day_run operations_period_night_run operations_period_mix_run PROPERTIES FIXTURES_SETUP ops_period_runs)
add_test(NAME operations_period_day_is_clocked_day COMMAND ${CMAKE_COMMAND} "-DPAIRS=integrated/ops_pday_dnl.csv;integrated/metrics_day_dnl.csv;integrated/ops_pday_lden.csv;integrated/metrics_day_lden.csv;integrated/ops_pday_cnel.csv;integrated/metrics_day_cnel.csv" -P tests/compare_pairs.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_period_night_matches_clock_dnl COMMAND golden_compare integrated/ops_pnight_dnl.csv integrated/metrics_2330_dnl.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_period_night_matches_clock_lden COMMAND golden_compare integrated/ops_pnight_lden.csv integrated/metrics_2330_lden.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_period_night_matches_clock_cnel COMMAND golden_compare integrated/ops_pnight_cnel.csv integrated/metrics_2330_cnel.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_period_day_is_clocked_day operations_period_night_matches_clock_dnl operations_period_night_matches_clock_lden operations_period_night_matches_clock_cnel PROPERTIES FIXTURES_REQUIRED "ops_period_runs;metrics_runs")
add_test(NAME operations_period_lden_recombines COMMAND golden_compare integrated/ops_pmix_lden.csv integrated/ops_pday_lden.csv 0 11.807644454215405 2 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_period_dnl_counts_evening_as_day COMMAND golden_compare integrated/ops_pmix_dnl.csv integrated/ops_pday_dnl.csv 0 11.139433523068368 2 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_period_cnel_recombines COMMAND golden_compare integrated/ops_pmix_cnel.csv integrated/ops_pday_cnel.csv 0 11.760912590556813 2 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_period_lden_recombines operations_period_dnl_counts_evening_as_day operations_period_cnel_recombines PROPERTIES FIXTURES_REQUIRED ops_period_runs)
add_test(NAME operations_period_clockless_runs COMMAND niceops --operations tests/data/operations/noclock_periods.csv -d databases/S-76D_M3.nod -v databases/S-76D.json --origin 39.82914323 -74.40996141 --frame_bearing 210 --ground_elevation 54.2651 --gridstep 1640 -D integrated/ops_noclock_periods_dnl.csv -o integrated/ops_noclock_periods.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_period_clockless_runs PROPERTIES FAIL_REGULAR_EXPRESSION "needs the time of day;Skipping" FIXTURES_SETUP ops_noclock_period_run)
# Its DNL is the night's, not an empty one: the same track dated 01:00 by --start_utc.
add_test(NAME operations_period_clockless_night_run COMMAND niceops -t tests/data/demo_geodetic_noclock.csv -d databases/S-76D_M3.nod -v databases/S-76D.json --origin 39.82914323 -74.40996141 --frame_bearing 210 --ground_elevation 54.2651 --gridstep 1640 --start_utc 3600 -D integrated/ops_noclock_night_dnl.csv -o integrated/ops_noclock_night.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_period_clockless_night_run PROPERTIES FIXTURES_SETUP ops_noclock_period_run)
add_test(NAME operations_period_clockless_is_the_night COMMAND golden_compare integrated/ops_noclock_periods_dnl.csv integrated/ops_noclock_night_dnl.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_period_clockless_is_the_night PROPERTIES FIXTURES_REQUIRED ops_noclock_period_run)
# A fleet: rows A, B, A on two databases, against -t A -t B -t A with per-track -d and -v.  Each
# database loads once, and each row's SEL file is named after it.
add_test(NAME operations_fleet_run COMMAND sh -c "rm -f integrated/ops_fleet.log && $<TARGET_FILE:niceops> --no-spectral --no-ground_effect --gridstep 2000 -L info -l integrated/ops_fleet.log --operations tests/data/operations/fleet.csv -D integrated/ops_fleet_dnl.csv -o integrated/ops_fleet/sel.csv" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_fleet_repeated_t COMMAND niceops --no-spectral --no-ground_effect --gridstep 2000 -L error -t tracks/S-76D/197421AC.csv -t tracks/S-76D/197422AC.csv -t tracks/S-76D/197421AC.csv -d databases/S-76D_M3.nod -d databases/B407_fixed_load_factor.nod -d databases/S-76D_M3.nod -v databases/S-76D.json -v databases/B407.json -v databases/S-76D.json -D integrated/ops_fleet_rep_dnl.csv -o integrated/ops_fleet_rep_a1.csv -o integrated/ops_fleet_rep_b.csv -o integrated/ops_fleet_rep_a2.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_fleet_run operations_fleet_repeated_t PROPERTIES FIXTURES_SETUP ops_fleet_runs)
add_test(NAME operations_fleet_mix COMMAND ${CMAKE_COMMAND} "-DPAIRS=integrated/ops_fleet_dnl.csv;integrated/ops_fleet_rep_dnl.csv;integrated/ops_fleet/a1_LAE.csv;integrated/ops_fleet_rep_a1.csv;integrated/ops_fleet/b_LAE.csv;integrated/ops_fleet_rep_b.csv;integrated/ops_fleet/a2_LAE.csv;integrated/ops_fleet_rep_a2.csv" -P tests/compare_pairs.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_fleet_loads_each_once COMMAND ${CMAKE_COMMAND} -E cat integrated/ops_fleet.log WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_fleet_mix operations_fleet_loads_each_once PROPERTIES FIXTURES_REQUIRED ops_fleet_runs)
set_tests_properties(operations_fleet_loads_each_once PROPERTIES PASS_REGULAR_EXPRESSION "Loading database databases/B407_fixed_load_factor\\.nod for tracks/S-76D/197422AC\\.csv" FAIL_REGULAR_EXPRESSION "Opening netCDF file: databases/S-76D_M3\\.nod.*Opening netCDF file: databases/S-76D_M3\\.nod;Opening netCDF file: databases/B407_fixed_load_factor\\.nod.*Opening netCDF file: databases/B407_fixed_load_factor\\.nod")
# The rows' resolved paths are read before the run starts, and still reach -l.
add_test(NAME operations_fleet_logs_its_paths COMMAND ${CMAKE_COMMAND} -E cat integrated/ops_fleet.log WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_fleet_logs_its_paths PROPERTIES FIXTURES_REQUIRED ops_fleet_runs PASS_REGULAR_EXPRESSION "\\[info\\] Operation a1 \\(tests/data/operations/fleet\\.csv line 2\\): track tracks/S-76D/197421AC\\.csv")
# The manifest records the operations file beside the tracks it names.
add_test(NAME operations_manifest_run COMMAND niceops --no-spectral --no-ground_effect --gridstep 2000 -L error -d databases/S-76D_M3.nod -v databases/S-76D.json --operations tests/data/operations/count3.csv -o integrated/ops_manifest.csv --manifest integrated/ops_manifest.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_manifest_run PROPERTIES FIXTURES_SETUP ops_manifest_run)
add_test(NAME operations_manifest_records_the_file COMMAND ${CMAKE_COMMAND} -E cat integrated/ops_manifest.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_manifest_records_the_file PROPERTIES FIXTURES_REQUIRED ops_manifest_run PASS_REGULAR_EXPRESSION "\"role\": \"operations\",[^}]*\"path\": \"tests/data/operations/count3\\.csv\",[^}]*\"sha256\": \"[0-9a-f]+\".*\"role\": \"trajectory\",[^}]*\"path\": \"tracks/S-76D/197421AC\\.csv\"")
# Refusals, by their messages.
set(OPS_DIR tests/data/operations)
add_test(NAME operations_refuses_t COMMAND niceops --operations ${OPS_DIR}/count1.csv -t tracks/S-76D/197421AC.csv -d databases/S-76D_M3.nod -v databases/S-76D.json -o integrated/ops_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_refuses_t PROPERTIES PASS_REGULAR_EXPRESSION "give them in its track column, not with -t")
add_test(NAME operations_refuses_two_d COMMAND niceops --operations ${OPS_DIR}/count1.csv -d databases/S-76D_M3.nod -d databases/S-76D_fixed_load_factor.nod -v databases/S-76D.json -o integrated/ops_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_refuses_two_d PROPERTIES PASS_REGULAR_EXPRESSION "With --operations, give per-operation -d in the file's database column")
add_test(NAME operations_needs_a_database COMMAND niceops --operations ${OPS_DIR}/no_database.csv -v databases/S-76D.json -o integrated/ops_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_needs_a_database PROPERTIES PASS_REGULAR_EXPRESSION "operation 197421AC \\(tests/data/operations/no_database\\.csv line 2\\) names no database")
add_test(NAME operations_missing_track COMMAND niceops --operations ${OPS_DIR}/missing_track.csv -d databases/S-76D_M3.nod -v databases/S-76D.json -o integrated/ops_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_missing_track PROPERTIES PASS_REGULAR_EXPRESSION "line 2 \\(no_such_track\\): the track file tests/data/operations/no_such_track\\.csv does not exist")
# A refusal of the operations file reaches -l, as every refusal does.
add_test(NAME operations_refusal_run COMMAND sh -c "rm -f integrated/ops_refusal.log; $<TARGET_FILE:niceops> --operations ${OPS_DIR}/missing_track.csv -d databases/S-76D_M3.nod -v databases/S-76D.json -o integrated/ops_refused.csv -l integrated/ops_refusal.log; test -f integrated/ops_refusal.log" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_refusal_run PROPERTIES FIXTURES_SETUP ops_refusal_run)
add_test(NAME operations_refusal_is_logged COMMAND ${CMAKE_COMMAND} -E cat integrated/ops_refusal.log WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_refusal_is_logged PROPERTIES FIXTURES_REQUIRED ops_refusal_run PASS_REGULAR_EXPRESSION "\\[error\\] [^\n]*the track file tests/data/operations/no_such_track\\.csv does not exist")
add_test(NAME operations_duplicate_names COMMAND niceops --operations ${OPS_DIR}/duplicate_names.csv -d databases/S-76D_M3.nod -v databases/S-76D.json -o integrated/ops_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_duplicate_names PROPERTIES PASS_REGULAR_EXPRESSION "the rows on lines 2 and 3 share the name 197421AC; add a name column")
add_test(NAME operations_bare_wind_speed COMMAND niceops ${OPS_RUN} --operations ${OPS_DIR}/bare_wind.csv -o integrated/ops_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_bare_wind_speed PROPERTIES PASS_REGULAR_EXPRESSION "--wind: name the speed's unit")
add_test(NAME operations_empty_air_cell COMMAND niceops -d databases/S-76D_M3.nod -v databases/S-76D.json --operations ${OPS_DIR}/empty_air.csv -o integrated/ops_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_empty_air_cell PROPERTIES PASS_REGULAR_EXPRESSION "the humidity column is empty\\. Fill it in every row, or give --humidity once")
add_test(NAME operations_offset_unit_mismatch COMMAND niceops --operations ${OPS_DIR}/offset_meters.csv -d databases/S-76D_M3.nod -v databases/S-76D.json -o integrated/ops_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_offset_unit_mismatch PROPERTIES PASS_REGULAR_EXPRESSION "the offset column declares m, and this run reads lengths in ft")
add_test(NAME operations_count_with_periods COMMAND niceops --operations ${OPS_DIR}/count_and_day.csv -d databases/S-76D_M3.nod -v databases/S-76D.json -o integrated/ops_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_count_with_periods PROPERTIES PASS_REGULAR_EXPRESSION "count_and_day\\.csv: give each row's operations as count or as day, evening and night counts, not both")
add_test(NAME operations_local_time_with_periods COMMAND niceops --operations ${OPS_DIR}/local_time_periods.csv -d databases/S-76D_M3.nod -v databases/S-76D.json -o integrated/ops_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_local_time_with_periods PROPERTIES PASS_REGULAR_EXPRESSION "local_time_periods\\.csv: local_time dates an operation[^\n]*give one or the other")
add_test(NAME operations_laeq_with_periods COMMAND niceops --operations ${OPS_DIR}/period_day.csv -d databases/S-76D_M3.nod -v databases/S-76D.json --laeq_file integrated/ops_refused_laeq.csv --laeq_window 7 19 -o integrated/ops_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_laeq_with_periods PROPERTIES PASS_REGULAR_EXPRESSION "--laeq_file is the level over a clock window, and tests/data/operations/period_day\\.csv gives day, evening and night counts")
# Lateral subtracks.  An offset of 0 (here under a unit token that agrees with
# --units) changes nothing, on an AC track and on a geodetic one.
add_test(NAME operations_offset_zero_run COMMAND niceops --no-spectral --no-ground_effect -d databases/S-76D_M3.nod -v databases/S-76D.json --operations ${OPS_DIR}/offset_zero.csv ${DEMO_COARSE} -o integrated/ops_offset_zero.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_offset_zero_run PROPERTIES FIXTURES_SETUP ops_offset_zero_run)
add_test(NAME operations_offset_zero_is_bit_identical COMMAND ${CMAKE_COMMAND} -E compare_files integrated/ops_offset_zero.csv integrated/demo_coarse.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_offset_zero_is_bit_identical PROPERTIES FIXTURES_REQUIRED "ops_offset_zero_run;demo_coarse_run")
set(OPS_GEODETIC -d databases/S-76D_M3.nod -v databases/S-76D.json --no-spectral --no-ground_effect --origin 39.82914323 -74.40996141 --frame_bearing 210 --ground_elevation 54.2651 --gridstep 1640 -L warn)
add_test(NAME operations_offset_zero_geodetic_run COMMAND niceops ${OPS_GEODETIC} --operations ${OPS_DIR}/offset_geodetic_zero.csv -o integrated/ops_offset_zero_geodetic.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_offset_zero_geodetic_t COMMAND niceops ${OPS_GEODETIC} -t tests/data/demo_geodetic_noclock.csv -o integrated/ops_offset_zero_geodetic_t.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_offset_zero_geodetic_run operations_offset_zero_geodetic_t PROPERTIES FIXTURES_SETUP ops_offset_zero_geodetic_runs)
add_test(NAME operations_offset_zero_geodetic_is_bit_identical COMMAND ${CMAKE_COMMAND} -E compare_files integrated/ops_offset_zero_geodetic.csv integrated/ops_offset_zero_geodetic_t.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_offset_zero_geodetic_is_bit_identical PROPERTIES FIXTURES_REQUIRED ops_offset_zero_geodetic_runs)
# A straight eastbound geodetic track offset 500 ft left (and right) is the same
# track flown 500 ft north (and south), on a grid asymmetric in y so that a sign
# error shows.  The hidden Catch case [.offset_tracks] writes the three tracks.
add_test(NAME operations_offset_tracks COMMAND niceops_unit_tests "[.offset_tracks]" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_offset_tracks PROPERTIES FIXTURES_SETUP offset_tracks)
set(OPS_MIRROR -d databases/S-76D_M3.nod -v databases/S-76D.json --no-spectral --no-ground_effect --origin 40 -75 --grid -6000 6000 -2000 4000 --gridstep 500 -L warn)
add_test(NAME operations_offset_left_run COMMAND niceops ${OPS_MIRROR} --operations ${OPS_DIR}/offset_left.csv -o integrated/ops_offset_left.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_offset_left_flown COMMAND niceops ${OPS_MIRROR} -t integrated/offset_left.csv -o integrated/ops_offset_left_flown.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_offset_right_run COMMAND niceops ${OPS_MIRROR} --operations ${OPS_DIR}/offset_right.csv -o integrated/ops_offset_right.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_offset_right_flown COMMAND niceops ${OPS_MIRROR} -t integrated/offset_right.csv -o integrated/ops_offset_right_flown.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_offset_left_run operations_offset_left_flown operations_offset_right_run operations_offset_right_flown PROPERTIES FIXTURES_REQUIRED offset_tracks FIXTURES_SETUP ops_offset_mirror_runs)
add_test(NAME operations_offset_mirror_left COMMAND golden_compare integrated/ops_offset_left.csv integrated/ops_offset_left_flown.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_offset_mirror_right COMMAND golden_compare integrated/ops_offset_right.csv integrated/ops_offset_right_flown.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_offset_mirror_left operations_offset_mirror_right PROPERTIES FIXTURES_REQUIRED ops_offset_mirror_runs)
# An offset adds no gain, so the peak search's bound needs nothing new; this
# pins that a subtrack through the spectral peak search equals a full march.
add_test(NAME operations_offset_windowed COMMAND niceops ${OPS_RUN} --operations ${OPS_DIR}/offset.csv -o integrated/ops_offset_windowed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME operations_offset_marched COMMAND niceops ${OPS_RUN} --operations ${OPS_DIR}/offset.csv -D integrated/ops_offset_dnl.csv -o integrated/ops_offset_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_offset_windowed operations_offset_marched PROPERTIES FIXTURES_SETUP ops_offset_window_runs)
add_test(NAME operations_windowed_is_exact_with_offset COMMAND golden_compare integrated/ops_offset_windowed.csv integrated/ops_offset_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_windowed_is_exact_with_offset PROPERTIES FIXTURES_REQUIRED ops_offset_window_runs)
# And it moves the footprint: the subtrack is not the centerline.
add_test(NAME operations_offset_moves_the_footprint COMMAND golden_compare integrated/ops_offset_windowed.csv integrated/ops_t.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_offset_moves_the_footprint PROPERTIES FIXTURES_REQUIRED "ops_offset_window_runs;ops_totals_runs" WILL_FAIL TRUE)
add_test(NAME operations_offset_on_a_static_hold COMMAND niceops ${OPS_RUN} --operations ${OPS_DIR}/offset_static.csv -o integrated/ops_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_offset_on_a_static_hold PROPERTIES PASS_REGULAR_EXPRESSION "static_hover\\.csv is a static operation, which has no direction of travel for a lateral offset")
add_test(NAME operations_offset_folds COMMAND niceops ${OPS_RUN} --operations ${OPS_DIR}/offset_fold.csv -o integrated/ops_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_offset_folds PROPERTIES PASS_REGULAR_EXPRESSION "a lateral offset of 1e\\+06 ft folds the path at t = [^\n]*must stay under 90% of its radius")
add_test(NAME operations_offset_needs_travel COMMAND niceops ${OPS_RUN} --operations ${OPS_DIR}/offset_hover.csv -o integrated/ops_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_offset_needs_travel PROPERTIES PASS_REGULAR_EXPRESSION "hover_geodetic\\.csv never moves faster than 3 kt over the ground")
add_test(NAME operations_provenance_two_databases COMMAND niceops --operations ${OPS_DIR}/fleet.csv -o integrated/ops_refused.csv --provenance_file integrated/ops_refused_provenance.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_provenance_two_databases PROPERTIES PASS_REGULAR_EXPRESSION "run it once per database")
add_test(NAME operations_clockless_needs_local_time COMMAND niceops --operations ${OPS_DIR}/noclock.csv -d databases/S-76D_M3.nod -v databases/S-76D.json --origin 39.82914323 -74.40996141 --frame_bearing 210 --ground_elevation 54.2651 --gridstep 1640 -D integrated/ops_noclock_dnl.csv -o integrated/ops_noclock.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(operations_clockless_needs_local_time PROPERTIES PASS_REGULAR_EXPRESSION "Skipping operation demo_geodetic_noclock \\(tests/data/operations/noclock\\.csv line 2\\): -D needs the time of day[^\n]*Give --start_utc, or local_time in the operations file")
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
# The ground-model tables' common build: the demo database's air and ground, source heights
# spanning the demo track's (40-1130 ft), one receiver row at 1.2 m, 10 ft in range.
set(GROUND_TABLE_BUILD --database databases/B407_fixed_load_factor.nod --sound_speed 335.636 --flow_resistance 225
    --source_heights 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200 210 220 230 240 250 260 270 280 290 300 310 320 330 340 350 360 370 380 390 400 410 420 430 440 450 460 470 480 490 500 510 520 530 540 550 560 570 580 590 600 610 620 630 640 650 660 670 680 690 700 710 720 730 740 750 760 770 780 790 800 810 820 830 840 850 860 870 880 890 900 910 920 930 940 950 960 970 980 990 1000 1010 1020 1030 1040 1050 1060 1070 1080 1090 1100 1110 1120 1130
    --receiver_heights 3.93701 --max_range 15000 --step 10)
# --propagation_table, end to end: a table built by niceops_propagation_table
# from the straight-ray ground model must give the footprint --ground_effect
# gives, to within the table's interpolation (10 ft in height and range: 0.19 dB
# at worst on this grid, 0.02 at 5 ft).  It checks the generator, the file, the lookup and the
# footprint's use of it together.
add_test(NAME propagation_table_build COMMAND niceops_propagation_table --method ground_model ${GROUND_TABLE_BUILD} -o tables/ground_model_test.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_run COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --propagation_table tables/ground_model_test.nc --gridstep 492 -L warn -o integrated/propagation_table.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_reference COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --ground_effect --exact_ground --gridstep 492 -L warn -o integrated/propagation_table_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_build PROPERTIES FIXTURES_SETUP propagation_table_file)
set_tests_properties(propagation_table_run PROPERTIES FIXTURES_REQUIRED propagation_table_file FIXTURES_SETUP propagation_table_runs)
set_tests_properties(propagation_table_reference PROPERTIES FIXTURES_SETUP propagation_table_runs)
# The peak search with a table: its bound takes the table's largest gain per
# band, so a windowed run must equal a marched one (-D forces the march) to the bit.
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
# A table pair takes its turbulence and roughness from the table.  A ground-model table built with
# them reproduces --ground_effect with the same (0.1 dB at worst here, where turbulence and roughness
# move cells by over 1 dB), Chessell and HARMONOISE alike, without the mismatch warning; a run asking
# for turbulence its table lacks is warned that only the uncovered pairs get it.
add_test(NAME propagation_table_turbulent_build COMMAND niceops_propagation_table --method ground_model ${GROUND_TABLE_BUILD} --turbulence 4e-4 --roughness 0.3 -o tables/ground_model_turbulent_test.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_harmonoise_build COMMAND niceops_propagation_table --method ground_model ${GROUND_TABLE_BUILD} --gamma_t 3e-5 -o tables/ground_model_harmonoise_test.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_turbulent_build propagation_table_harmonoise_build PROPERTIES FIXTURES_SETUP propagation_table_turbulent_files)
add_test(NAME propagation_table_turbulent_run COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --propagation_table tables/ground_model_turbulent_test.nc --turbulence 4e-4 --roughness 0.3 --gridstep 492 -L warn -o integrated/propagation_table_turbulent.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_turbulent_reference COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --ground_effect --exact_ground --turbulence 4e-4 --roughness 0.3 --gridstep 492 -L warn -o integrated/propagation_table_turbulent_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_harmonoise_run COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --propagation_table tables/ground_model_harmonoise_test.nc --turbulence_model harmonoise --gamma_t 3e-5 --gridstep 492 -L warn -o integrated/propagation_table_harmonoise.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_harmonoise_reference COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --ground_effect --exact_ground --turbulence_model harmonoise --gamma_t 3e-5 --gridstep 492 -L warn -o integrated/propagation_table_harmonoise_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_turbulent_run propagation_table_harmonoise_run PROPERTIES FIXTURES_REQUIRED propagation_table_turbulent_files FAIL_REGULAR_EXPRESSION "was built with")
set_tests_properties(propagation_table_turbulent_run propagation_table_turbulent_reference propagation_table_harmonoise_run propagation_table_harmonoise_reference PROPERTIES FIXTURES_SETUP propagation_table_turbulent_runs)
add_test(NAME propagation_table_carries_turbulence COMMAND golden_compare integrated/propagation_table_turbulent.csv integrated/propagation_table_turbulent_reference.csv 0.1 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_carries_harmonoise COMMAND golden_compare integrated/propagation_table_harmonoise.csv integrated/propagation_table_harmonoise_reference.csv 0.1 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME propagation_table_turbulence_matters COMMAND golden_compare integrated/propagation_table_turbulent_reference.csv integrated/propagation_table_reference.csv 0.5 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_carries_turbulence propagation_table_carries_harmonoise propagation_table_turbulence_matters PROPERTIES FIXTURES_REQUIRED "propagation_table_turbulent_runs;propagation_table_runs")
set_tests_properties(propagation_table_turbulence_matters PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value")
add_test(NAME propagation_table_warns_turbulence_mismatch COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --propagation_table tables/ground_model_test.nc --turbulence 4e-4 --gridstep 4000 -L warn -o integrated/propagation_table_mismatch.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_warns_turbulence_mismatch PROPERTIES FIXTURES_REQUIRED propagation_table_file PASS_REGULAR_EXPRESSION "was built with still air, roughness 0 m, and its pairs take that; the run's Chessell turbulence 0\\.0004 rad s/sqrt\\(m\\), roughness 0 m reaches only the pairs it does not cover")
add_test(NAME propagation_table_warns_still_run_on_turbulent_table COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --propagation_table tables/ground_model_turbulent_test.nc --gridstep 4000 -L warn -o integrated/propagation_table_mismatch_still.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_warns_still_run_on_turbulent_table PROPERTIES FIXTURES_REQUIRED propagation_table_turbulent_files PASS_REGULAR_EXPRESSION "the table was built with Chessell turbulence 0\\.0004 rad s/sqrt\\(m\\), roughness 0\\.09144 m, and its pairs take that; the run's still air, roughness 0 m")
add_test(NAME propagation_table_turbulence_refusals COMMAND niceops_propagation_table --method ground_model --bands 100 --source_heights 100 --receiver_heights 4 --max_range 1000 --step 100 --turbulence 4e-4 --gamma_t 1e-5 -o integrated/never.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(propagation_table_turbulence_refusals PROPERTIES PASS_REGULAR_EXPRESSION "give the ground model one coherence")
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
# --urban_table (urban.hpp), end to end: the demo track on the B407 database, at seven microphones
# (tests/data/urban_observers.csv), poles and flush plates, in ideal canyons built by
# niceops_urban_table (tests/data/urban_canyons.csv).  A table of open microphones -- 0 dB and no
# facades -- is the run without one, byte for byte.  The peak search's bound takes the table's
# largest gain and each first-order image in its own direction, so a windowed run equals its march
# (-D forces it), for the canyons, on plates, and for a table with +25 dB at 3 to 20 degrees of
# elevation, which kills a bound that leaves the table's gain out.  The image term of the bound has
# no CTest that fails without it: on this database the table's own largest gain and the bound's
# slack already cover the images (docs/internals/testing.md).
set(URBAN_RUN niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --frame_bearing 210 -b tests/data/urban_observers.csv -L warn)
# The runs over the test district (--buildings tests/data/urban_buildings.csv) were two thirds of the
# suite's time on the full demo track.  The states sampled against every state, which needs whole
# passes, take the track at every 8th state (tests/data/urban_track_every8.csv): each microphone keeps
# its pass, the SELs stay within 0.03 dB of the full track's, and the sampling still evaluates under
# half the pairs.  Every other district run checks an
# option against the default and takes 4 s of the approach (URBAN_DISTRICT_SHORT), a tenth again, as do
# the far building's degenerate cases below, which hold state by state.
set(URBAN_DISTRICT niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tests/data/urban_track_every8.csv --frame_bearing 210 -b tests/data/urban_observers.csv -L warn)
set(URBAN_DISTRICT_SHORT niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tests/data/urban_track_short.csv --frame_bearing 210 -b tests/data/urban_observers.csv -L warn)
set(URBAN_SHORT_BARE niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tests/data/urban_track_short.csv --frame_bearing 210 -b tests/data/urban_shadow_observers.csv -L warn)
set(URBAN_PLATE --flow_resistance 225 --ground_effect --plate_table plates/gr1425_db225.nc)
add_test(NAME urban_table_build COMMAND niceops_urban_table --canyons tests/data/urban_canyons.csv --database databases/B407_fixed_load_factor.nod --frame_bearing 210 -o tables/urban_canyon_test.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_open_build COMMAND niceops_urban_table --canyons tests/data/urban_open.csv --database databases/B407_fixed_load_factor.nod --frame_bearing 210 -o tables/urban_open_test.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_build urban_table_open_build PROPERTIES FIXTURES_SETUP urban_table_files)
add_test(NAME urban_table_bump_build COMMAND niceops_unit_tests "[.urban_bump_table]" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_bump_build PROPERTIES FIXTURES_REQUIRED urban_table_files FIXTURES_SETUP urban_table_bump_file)
add_test(NAME urban_table_reference COMMAND ${URBAN_RUN} --ground_effect -o integrated/urban_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_plate_reference COMMAND ${URBAN_RUN} ${URBAN_PLATE} -o integrated/urban_plate_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_reference urban_table_plate_reference PROPERTIES FIXTURES_SETUP urban_table_runs)
add_test(NAME urban_table_open COMMAND ${URBAN_RUN} --ground_effect --urban_table tables/urban_open_test.nc -o integrated/urban_open.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_canyon COMMAND ${URBAN_RUN} --ground_effect --urban_table tables/urban_canyon_test.nc -o integrated/urban_canyon.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_canyon_marched COMMAND ${URBAN_RUN} --ground_effect --urban_table tables/urban_canyon_test.nc -D integrated/urban_canyon_dnl.csv -o integrated/urban_canyon_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_plate_open COMMAND ${URBAN_RUN} ${URBAN_PLATE} --urban_table tables/urban_open_test.nc -o integrated/urban_plate_open.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_plate COMMAND ${URBAN_RUN} ${URBAN_PLATE} --urban_table tables/urban_canyon_test.nc -o integrated/urban_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_plate_marched COMMAND ${URBAN_RUN} ${URBAN_PLATE} --urban_table tables/urban_canyon_test.nc -D integrated/urban_plate_dnl.csv -o integrated/urban_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_open urban_table_canyon urban_table_canyon_marched urban_table_plate_open urban_table_plate urban_table_plate_marched PROPERTIES FIXTURES_REQUIRED urban_table_files FIXTURES_SETUP urban_table_runs)
add_test(NAME urban_table_bump COMMAND ${URBAN_RUN} --ground_effect --urban_table tables/urban_bump_test.nc -o integrated/urban_bump.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_bump_marched COMMAND ${URBAN_RUN} --ground_effect --urban_table tables/urban_bump_test.nc -D integrated/urban_bump_dnl.csv -o integrated/urban_bump_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_bump urban_table_bump_marched PROPERTIES FIXTURES_REQUIRED urban_table_bump_file FIXTURES_SETUP urban_table_runs)
add_test(NAME urban_table_open_is_open_ground COMMAND ${CMAKE_COMMAND} -E compare_files integrated/urban_open.csv integrated/urban_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_plate_open_is_open_ground COMMAND ${CMAKE_COMMAND} -E compare_files integrated/urban_plate_open.csv integrated/urban_plate_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_windowed_is_exact COMMAND golden_compare integrated/urban_canyon.csv integrated/urban_canyon_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_plate_windowed_is_exact COMMAND golden_compare integrated/urban_plate.csv integrated/urban_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_bump_windowed_is_exact COMMAND golden_compare integrated/urban_bump.csv integrated/urban_bump_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# The canyons move the footprint: by -3.7 to +2.0 dB at these seven microphones.
add_test(NAME urban_table_changes_the_footprint COMMAND golden_compare integrated/urban_canyon.csv integrated/urban_reference.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_open_is_open_ground urban_table_plate_open_is_open_ground urban_table_windowed_is_exact urban_table_plate_windowed_is_exact urban_table_bump_windowed_is_exact urban_table_changes_the_footprint PROPERTIES FIXTURES_REQUIRED urban_table_runs)
set_tests_properties(urban_table_changes_the_footprint PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value")
# The windowed/marched pairs compare a run with itself, so the canyon and plate runs are pinned
# (golden/urban_canyon.csv, golden/urban_plate.csv, made by the build that added --urban_table):
# an image path whose plate azimuth is not mirrored back moves the plate microphones.
add_test(NAME urban_table_canyon_golden COMMAND golden_compare integrated/urban_canyon.csv golden/urban_canyon.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_plate_golden COMMAND golden_compare integrated/urban_plate.csv golden/urban_plate.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_canyon_golden urban_table_plate_golden PROPERTIES FIXTURES_REQUIRED urban_table_runs)
# The wave band (--wave_below, src/canyon_wave): the flush canyons of the same file solved by the wave
# equation below 80 Hz and blended through 125 Hz (coarse, 10 points per wavelength, so it takes a
# second), the rest warned about and left geometric.  The windowed run equals its march, and the run
# is pinned (golden/urban_wave.csv).  Leaving the facades' reflections unscaled in the wave band
# moves the A-weighted SEL by under 1e-5 dB below 125 Hz, so the Catch2 case "the wave band's
# weights and blend" is what guards that.
add_test(NAME urban_table_wave_build COMMAND niceops_urban_table --canyons tests/data/urban_canyons.csv --database databases/B407_fixed_load_factor.nod --frame_bearing 210 --wave_below 80 --wave_ppw 10 -j 2 --no-gpu -o tables/urban_wave_test.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_wave_build PROPERTIES FIXTURES_SETUP urban_table_wave_file)
add_test(NAME urban_table_wave_warns_unequal_walls COMMAND niceops_urban_table --canyons tests/data/urban_canyons.csv --bands 50 63 80 100 --frame_bearing 210 --wave_below 63 --wave_ppw 8 -o tables/urban_wave_warn_test.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_wave_warns_unequal_walls PROPERTIES PASS_REGULAR_EXPRESSION "row 2: walls of unequal height \\(40 and 60 ft\\); the wave solver needs flush roofs")
add_test(NAME urban_table_wave COMMAND ${URBAN_RUN} --ground_effect --urban_table tables/urban_wave_test.nc -o integrated/urban_wave.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_wave_marched COMMAND ${URBAN_RUN} --ground_effect --urban_table tables/urban_wave_test.nc -D integrated/urban_wave_dnl.csv -o integrated/urban_wave_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_wave urban_table_wave_marched PROPERTIES FIXTURES_REQUIRED urban_table_wave_file FIXTURES_SETUP urban_table_wave_runs)
add_test(NAME urban_table_wave_windowed_is_exact COMMAND golden_compare integrated/urban_wave.csv integrated/urban_wave_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_table_wave_golden COMMAND golden_compare integrated/urban_wave.csv golden/urban_wave.csv 1e-5 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_wave_windowed_is_exact urban_table_wave_golden PROPERTIES FIXTURES_REQUIRED urban_table_wave_runs)
# --gpu marches the wave equation and nothing else, so without --wave_below it
# is refused before any work, whatever the build's backend.
add_test(NAME urban_table_gpu_without_wave_refused COMMAND niceops_urban_table --canyons tests/data/urban_canyons.csv --bands 50 63 80 --frame_bearing 210 --gpu -o integrated/never.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_gpu_without_wave_refused PROPERTIES PASS_REGULAR_EXPRESSION "there is none to march")
# The wave band marched on the GPU (--gpu, in float): the same SEL as the CPU's
# table within 0.001 dB.  A build without a GPU backend refuses the flag.
if(NICEOPS_GPU_BACKEND STREQUAL "OFF")
  add_test(NAME urban_table_wave_gpu_refused COMMAND niceops_urban_table --canyons tests/data/urban_canyons.csv --bands 50 63 80 --frame_bearing 210 --wave_below 63 --wave_ppw 8 --gpu -o integrated/never.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(urban_table_wave_gpu_refused PROPERTIES PASS_REGULAR_EXPRESSION "this build has no GPU backend")
else()
  add_test(NAME urban_table_wave_gpu_build COMMAND niceops_urban_table --canyons tests/data/urban_canyons.csv --database databases/B407_fixed_load_factor.nod --frame_bearing 210 --wave_below 80 --wave_ppw 10 -j 2 --gpu -o tables/urban_wave_gpu_test.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(urban_table_wave_gpu_build PROPERTIES FIXTURES_SETUP urban_table_wave_gpu_file)
  add_test(NAME urban_table_wave_gpu COMMAND ${URBAN_RUN} --ground_effect --urban_table tables/urban_wave_gpu_test.nc -o integrated/urban_wave_gpu.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(urban_table_wave_gpu PROPERTIES FIXTURES_REQUIRED urban_table_wave_gpu_file FIXTURES_SETUP urban_table_wave_gpu_run)
  add_test(NAME urban_table_wave_gpu_matches_cpu COMMAND golden_compare integrated/urban_wave_gpu.csv integrated/urban_wave.csv 0.001 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(urban_table_wave_gpu_matches_cpu PROPERTIES FIXTURES_REQUIRED "urban_table_wave_gpu_run;urban_table_wave_runs")
endif()
# LoD1 buildings (--buildings, src/building_paths): the same seven microphones among eleven buildings
# (tests/data/urban_buildings.csv) -- a street, an L's corner, an intersection, a wall beside a flush
# microphone, an open one, an alley and a building between a pole and the track -- as poles and as
# plates.  Each windowed run equals its march, and both are pinned (golden/urban_buildings.csv,
# golden/urban_buildings_plate.csv, made by the build that added --buildings).
add_test(NAME urban_buildings_build COMMAND niceops_urban_table --buildings tests/data/urban_buildings.csv --receivers tests/data/urban_observers.csv --database databases/B407_fixed_load_factor.nod --frame_bearing 210 -j 2 -o tables/urban_buildings_test.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings_build PROPERTIES FIXTURES_SETUP urban_buildings_file)
add_test(NAME urban_buildings COMMAND ${URBAN_RUN} --ground_effect --urban_table tables/urban_buildings_test.nc -o integrated/urban_buildings.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_buildings_marched COMMAND ${URBAN_RUN} --ground_effect --urban_table tables/urban_buildings_test.nc -D integrated/urban_buildings_dnl.csv -o integrated/urban_buildings_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_buildings_plate COMMAND ${URBAN_RUN} ${URBAN_PLATE} --urban_table tables/urban_buildings_test.nc -o integrated/urban_buildings_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_buildings_plate_marched COMMAND ${URBAN_RUN} ${URBAN_PLATE} --urban_table tables/urban_buildings_test.nc -D integrated/urban_buildings_plate_dnl.csv -o integrated/urban_buildings_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings urban_buildings_marched urban_buildings_plate urban_buildings_plate_marched PROPERTIES FIXTURES_REQUIRED urban_buildings_file FIXTURES_SETUP urban_buildings_runs)
add_test(NAME urban_buildings_windowed_is_exact COMMAND golden_compare integrated/urban_buildings.csv integrated/urban_buildings_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_buildings_plate_windowed_is_exact COMMAND golden_compare integrated/urban_buildings_plate.csv integrated/urban_buildings_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_buildings_golden COMMAND golden_compare integrated/urban_buildings.csv golden/urban_buildings.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_buildings_plate_golden COMMAND golden_compare integrated/urban_buildings_plate.csv golden/urban_buildings_plate.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings_windowed_is_exact urban_buildings_plate_windowed_is_exact urban_buildings_golden urban_buildings_plate_golden PROPERTIES FIXTURES_REQUIRED urban_buildings_runs)
# The same with walls scattering a fifth of what they reflect (--facade_scattering 0.2, src/diffuse_facades):
# the windowed run equals its march, which the diffuse field's share of the table's largest gain
# has to keep so, and the run is pinned (golden/urban_buildings_diffuse.csv).
add_test(NAME urban_buildings_diffuse_build COMMAND niceops_urban_table --buildings tests/data/urban_buildings.csv --receivers tests/data/urban_observers.csv --database databases/B407_fixed_load_factor.nod --frame_bearing 210 --facade_scattering 0.2 -j 2 -o tables/urban_buildings_diffuse_test.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings_diffuse_build PROPERTIES FIXTURES_SETUP urban_buildings_diffuse_file)
add_test(NAME urban_buildings_diffuse COMMAND ${URBAN_RUN} --ground_effect --urban_table tables/urban_buildings_diffuse_test.nc -o integrated/urban_buildings_diffuse.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_buildings_diffuse_marched COMMAND ${URBAN_RUN} --ground_effect --urban_table tables/urban_buildings_diffuse_test.nc -D integrated/urban_buildings_diffuse_dnl.csv -o integrated/urban_buildings_diffuse_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings_diffuse urban_buildings_diffuse_marched PROPERTIES FIXTURES_REQUIRED urban_buildings_diffuse_file FIXTURES_SETUP urban_buildings_diffuse_runs)
add_test(NAME urban_buildings_diffuse_windowed_is_exact COMMAND golden_compare integrated/urban_buildings_diffuse.csv integrated/urban_buildings_diffuse_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_buildings_diffuse_golden COMMAND golden_compare integrated/urban_buildings_diffuse.csv golden/urban_buildings_diffuse.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings_diffuse_windowed_is_exact urban_buildings_diffuse_golden PROPERTIES FIXTURES_REQUIRED urban_buildings_diffuse_runs)
# GeoJSON (src/geojson_buildings): the same district as WGS84 polygons (tests/data/urban_buildings.geojson,
# made from the CSV through the demo frame's inverse projection) gives the same run, to round-off; a
# GeoJSON file needs the run's --origin, and buildings without a height need --default_height.
add_test(NAME urban_buildings_geojson_build COMMAND niceops_urban_table --buildings tests/data/urban_buildings.geojson --origin 39.8104276 -74.4263086 --receivers tests/data/urban_observers.csv --database databases/B407_fixed_load_factor.nod --frame_bearing 210 -j 2 -o tables/urban_buildings_geojson_test.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings_geojson_build PROPERTIES FIXTURES_SETUP urban_buildings_geojson_file)
add_test(NAME urban_buildings_geojson COMMAND ${URBAN_RUN} --ground_effect --urban_table tables/urban_buildings_geojson_test.nc -o integrated/urban_buildings_geojson.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings_geojson PROPERTIES FIXTURES_REQUIRED urban_buildings_geojson_file FIXTURES_SETUP urban_buildings_geojson_runs)
add_test(NAME urban_buildings_geojson_is_the_csv COMMAND golden_compare integrated/urban_buildings_geojson.csv integrated/urban_buildings.csv 1e-7 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings_geojson_is_the_csv PROPERTIES FIXTURES_REQUIRED "urban_buildings_geojson_runs;urban_buildings_runs")
add_test(NAME urban_buildings_geojson_needs_origin COMMAND niceops_urban_table --buildings tests/data/urban_buildings.geojson --receivers tests/data/urban_observers.csv --bands 1000 -o integrated/never.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings_geojson_needs_origin PROPERTIES PASS_REGULAR_EXPRESSION "is GeoJSON, in longitude and latitude: give the run's --origin")
add_test(NAME urban_buildings_geojson_needs_heights COMMAND niceops_urban_table --buildings tests/data/urban_no_height.geojson --origin 39.8104276 -74.4263086 --receivers tests/data/urban_observers.csv --bands 1000 --frame_bearing 210 -o integrated/never.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings_geojson_needs_heights PROPERTIES PASS_REGULAR_EXPRESSION "1 buildings give no height \\(height, building:height or building:levels\\); give --default_height")
# The builder refuses a microphone inside a building, the wave band with buildings, and buildings
# without receivers.
# --buildings (urban_point.hpp): every microphone hears the source's own rays among the buildings.
# With one building 1000 km off (its reflections absorbed away, its edges far outside every ray's
# radius) and a floor that absorbs everything, the one ray left is the straight path, and the run is
# the run without a ground model: every state evaluated (--state_tolerance 0), at 1e-9; and with the
# states sampled, whose transfer is then nothing at every state, the direct ray's free-field level
# that the sampling adds it to is the ray's own, at 1e-6 dB.
add_test(NAME buildings_far_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_far_building.csv --floor_absorption 1 --phase_interval 90 --state_tolerance 0 --no-gpu -o integrated/buildings_far.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_far_sampled_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_far_building.csv --floor_absorption 1 --no-gpu -o integrated/buildings_far_sampled.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_far_reference COMMAND ${URBAN_DISTRICT_SHORT} --no-ground_effect -o integrated/buildings_far_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_far_run buildings_far_sampled_run buildings_far_reference PROPERTIES FIXTURES_SETUP buildings_far_runs)
add_test(NAME buildings_far_is_open_without_ground COMMAND golden_compare integrated/buildings_far.csv integrated/buildings_far_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_far_sampled_is_open COMMAND golden_compare integrated/buildings_far_sampled.csv integrated/buildings_far_reference.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_far_is_open_without_ground buildings_far_sampled_is_open PROPERTIES FIXTURES_REQUIRED buildings_far_runs)
# The demo approach over the test district: its states sampled (the default) against every state
# evaluated, within 0.01 dB of SEL at every microphone (0.001 measured; 0.003 on the full track).
add_test(NAME buildings_every_state_run COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --state_tolerance 0 --no-gpu -o integrated/buildings_every_state.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_sampled_run COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --no-gpu -o integrated/buildings_sampled.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_every_state_run buildings_sampled_run PROPERTIES FIXTURES_SETUP buildings_sampling_runs)
add_test(NAME buildings_states_sampled COMMAND golden_compare integrated/buildings_sampled.csv integrated/buildings_every_state.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_states_sampled PROPERTIES FIXTURES_REQUIRED buildings_sampling_runs)
# The option runs below take the district's default run over a short stretch of the track as their
# reference: what they check (an option that changes nothing, or something) needs no whole pass.
add_test(NAME buildings_short_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_buildings.csv --no-gpu -o integrated/buildings_short.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_short_run PROPERTIES FIXTURES_SETUP buildings_short_runs)
# --no-ground_effect is no clash (the run has no ground model among buildings); and a database
# without spectra is refused, its rays reading each band.  (An observer's z is its floor's height
# where it stands on no roof or wall: buildings_floor_at_observer_z below.)
add_test(NAME buildings_no_ground_effect_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_buildings.csv --no-ground_effect --no-gpu -o integrated/buildings_no_ground.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_no_ground_effect_run PROPERTIES FIXTURES_SETUP buildings_receiver_runs)
add_test(NAME buildings_no_ground_effect_is_the_default COMMAND golden_compare integrated/buildings_no_ground.csv integrated/buildings_short.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_no_ground_effect_is_the_default PROPERTIES FIXTURES_REQUIRED "buildings_receiver_runs;buildings_short_runs")
# --no_profile_pe moves nothing here: no line from the track to these microphones crosses a row of
# flat-roofed blocks with a street near grazing (the row's profile, src/profile_pe).
add_test(NAME buildings_no_profile_pe_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_buildings.csv --no_profile_pe --no-gpu -o integrated/buildings_no_profile_pe.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_no_profile_pe_run PROPERTIES FIXTURES_SETUP buildings_no_profile_runs)
add_test(NAME buildings_no_profile_pe_reaches_no_row COMMAND golden_compare integrated/buildings_no_profile_pe.csv integrated/buildings_short.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_no_profile_pe_reaches_no_row PROPERTIES FIXTURES_REQUIRED "buildings_no_profile_runs;buildings_short_runs")
# --path_radius is a fast mode: over the test district at 75 m (246 ft) rather than 150 the SEL moves by
# 0.0003 dB at most on the short track (0.009 on the full one, 0.014 there at 100 m); the default-radius
# run is its reference.  That alone would pass with the option ignored, so two more runs pin it: 150 m given in
# feet is the default run exactly (the radius is converted from --units: unconverted, 492 m changes the SELs), and 30 ft is not.
add_test(NAME buildings_path_radius_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_buildings.csv --path_radius 246 --no-gpu -o integrated/buildings_path_radius.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_path_radius_default_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_buildings.csv --path_radius 492.1259842519685 --no-gpu -o integrated/buildings_path_radius_default.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_path_radius_short_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_buildings.csv --path_radius 30 --no-gpu -o integrated/buildings_path_radius_short.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_path_radius_run buildings_path_radius_default_run buildings_path_radius_short_run PROPERTIES FIXTURES_SETUP buildings_path_radius_runs)
add_test(NAME buildings_path_radius_is_close COMMAND golden_compare integrated/buildings_path_radius.csv integrated/buildings_short.csv 0.03 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_path_radius_in_units COMMAND golden_compare integrated/buildings_path_radius_default.csv integrated/buildings_short.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_path_radius_takes_effect COMMAND golden_compare integrated/buildings_path_radius_short.csv integrated/buildings_short.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_path_radius_takes_effect PROPERTIES WILL_FAIL TRUE)
set_tests_properties(buildings_path_radius_is_close buildings_path_radius_in_units buildings_path_radius_takes_effect PROPERTIES FIXTURES_REQUIRED "buildings_path_radius_runs;buildings_short_runs")
add_test(NAME buildings_refuses_lone_path_radius COMMAND ${URBAN_RUN} --path_radius 246 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_lone_path_radius PROPERTIES PASS_REGULAR_EXPRESSION "--path_radius, --gpu and --no-gpu go with --buildings")
add_test(NAME buildings_refuses_zero_path_radius COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --path_radius 0 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_zero_path_radius PROPERTIES PASS_REGULAR_EXPRESSION "--path_radius is a positive length")
add_test(NAME buildings_refuses_no_spectra COMMAND niceops -d databases/S-76D_M3.nod -v databases/S-76D.json -t tracks/S-76D/197421AC.csv --frame_bearing 210 -b tests/data/urban_observers.csv -L warn --buildings tests/data/urban_buildings.csv -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_no_spectra PROPERTIES PASS_REGULAR_EXPRESSION "--buildings needs a database carrying per-band spectra")
add_test(NAME buildings_refuses_lone_state_tolerance COMMAND ${URBAN_RUN} --state_tolerance 0.1 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_lone_state_tolerance PROPERTIES PASS_REGULAR_EXPRESSION "--exact_shared_faces, --no_profile_pe, --tilted_step_split, --edge_source, --state_tolerance, --path_radius, --gpu and --no-gpu go with --buildings")
add_test(NAME buildings_refuses_lone_gpu COMMAND ${URBAN_RUN} --no-gpu -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_lone_gpu PROPERTIES PASS_REGULAR_EXPRESSION "--gpu and --no-gpu go with --buildings")
# Each microphone's floor at its observer's z where it stands on no roof or wall: observers 5 ft up
# hear what observers on the district's floor hear of the track 5 ft lower (a short stretch of the
# demo track and its copy lowered by 5 ft), to round-off; and a run with microphones on two floors
# goes through.
set(URBAN_SHORT niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 --frame_bearing 210 -L warn --buildings tests/data/urban_buildings.csv)
add_test(NAME buildings_floor_raised_run COMMAND ${URBAN_SHORT} -t tests/data/urban_track_short.csv -b tests/data/urban_observers_up5.csv -o integrated/urban_floor_raised.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_floor_raised_run PROPERTIES FIXTURES_SETUP urban_floor_raised)
add_test(NAME buildings_floor_lowered_run COMMAND ${URBAN_SHORT} -t tests/data/urban_track_short_lowered.csv -b tests/data/urban_observers.csv -o integrated/urban_floor_lowered.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_floor_lowered_run PROPERTIES FIXTURES_SETUP urban_floor_lowered)
add_test(NAME buildings_floor_at_observer_z COMMAND golden_compare integrated/urban_floor_raised.csv integrated/urban_floor_lowered.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_floor_at_observer_z PROPERTIES FIXTURES_REQUIRED "urban_floor_raised;urban_floor_lowered")
add_test(NAME buildings_two_floors_run COMMAND ${URBAN_SHORT} -t tests/data/urban_track_short.csv -b tests/data/urban_observers_raised.csv -o integrated/urban_two_floors.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# Pitched roofs (--roofs): the demo district's two terraces given gables runs through, and --roofs
# alone is refused.
add_test(NAME buildings_pitched_runs COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_buildings.csv --roofs tests/data/urban_roofs.csv -o integrated/urban_pitched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# A microphone flush on a slope of a gable (building 1's, 15 ft up the side from its eave) takes the slope's
# reflection and runs; with a plate, which lies level, it is refused.
set(URBAN_SLOPE niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tests/data/urban_track_short.csv --frame_bearing 210 -b tests/data/urban_on_slope.csv -L warn --buildings tests/data/urban_buildings.csv --roofs tests/data/urban_roofs.csv)
add_test(NAME buildings_slope_microphone_runs COMMAND ${URBAN_SLOPE} -o integrated/urban_on_slope.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_refuses_plate_on_slope COMMAND ${URBAN_SLOPE} ${URBAN_PLATE} -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_plate_on_slope PROPERTIES PASS_REGULAR_EXPRESSION "observer 1's microphone is flush on a wall or a sloping roof")
add_test(NAME buildings_refuses_lone_roofs COMMAND ${URBAN_RUN} --roofs tests/data/urban_roofs.csv -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_lone_roofs PROPERTIES PASS_REGULAR_EXPRESSION "--roofs goes with --buildings")
add_test(NAME buildings_refuses_negative_state_tolerance COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --state_tolerance -1 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_negative_state_tolerance PROPERTIES PASS_REGULAR_EXPRESSION "--state_tolerance is a non-negative number of dB")
# Facades that scatter (--facade_scattering, diffuse_facades.hpp): each microphone's diffuse field, lit
# by the source's specular paths with the sphere read along each one's launch, adds to its rays.  The
# states sampled against every state evaluated, within 0.01 dB of SEL (0.004 measured), and the run
# pinned (golden/buildings_scattering.csv, made by the build that added it).  A building outside every
# microphone's --diffuse_radius keeps its whole rho, so with none in reach the run is the specular one.
# A spectrum (--facade_scattering_spectrum) must reach every band the metric sums.  The library test
# `a run's diffuse field is niceops_point_field's for an omnidirectional source` holds the field to the
# tool's, on the CPU and the GPU.
add_test(NAME buildings_scattering_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_buildings.csv --facade_scattering 0.3 --patch_size 20 --no-gpu -o integrated/buildings_scattering.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_scattering_out_of_reach_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_buildings.csv --facade_scattering 0.5 --diffuse_radius 0.001 --no-gpu -o integrated/buildings_scattering_out_of_reach.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_scattering_spectrum_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_buildings.csv --facade_scattering_spectrum tests/data/scattering_spectrum.csv --band_limits 100 2000 --patch_size 20 --no-gpu -o integrated/buildings_scattering_spectrum.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_scattering_run buildings_scattering_out_of_reach_run buildings_scattering_spectrum_run PROPERTIES FIXTURES_SETUP buildings_scattering_runs)
add_test(NAME buildings_scattering_golden COMMAND golden_compare integrated/buildings_scattering.csv golden/buildings_scattering.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_scattering_out_of_reach_is_specular COMMAND golden_compare integrated/buildings_scattering_out_of_reach.csv integrated/buildings_short.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_scattering_golden buildings_scattering_out_of_reach_is_specular PROPERTIES FIXTURES_REQUIRED "buildings_scattering_runs;buildings_short_runs")
add_test(NAME buildings_refuses_scattering_out_of_range COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --facade_scattering 1.5 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_scattering_out_of_range PROPERTIES PASS_REGULAR_EXPRESSION "--facade_scattering is a share within")
add_test(NAME buildings_refuses_short_scattering_spectrum COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --facade_scattering_spectrum tests/data/scattering_spectrum.csv -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_short_scattering_spectrum PROPERTIES PASS_REGULAR_EXPRESSION "covers 100 to 2000 Hz, not the band at")
add_test(NAME buildings_refuses_lone_scattering COMMAND ${URBAN_RUN} --facade_scattering 0.3 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_lone_scattering PROPERTIES PASS_REGULAR_EXPRESSION "go with --buildings")
# Turbulence scattering into the buildings' shadows (--scattering under --buildings,
# turbulence_scattering.hpp): on the demo approach every microphone sees the helicopter at some
# point, and strong turbulence moves no SEL by 0.01 dB (0.007 measured); without --buildings or
# --refraction it is still refused.  The library tests hold the fit, the screen and the run.
add_test(NAME buildings_turbulence_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_buildings.csv --scattering 1e-5 --no-gpu -o integrated/buildings_turbulence.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_turbulence_run PROPERTIES FIXTURES_SETUP buildings_turbulence_runs)
add_test(NAME buildings_turbulence_is_small_on_the_demo COMMAND golden_compare integrated/buildings_turbulence.csv integrated/buildings_short.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_turbulence_is_small_on_the_demo PROPERTIES FIXTURES_REQUIRED "buildings_turbulence_runs;buildings_short_runs")
add_test(NAME scattering_refuses_open_ground_without_refraction COMMAND ${URBAN_RUN} --scattering 1e-6 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(scattering_refuses_open_ground_without_refraction PROPERTIES PASS_REGULAR_EXPRESSION "give --refraction and --ground_effect, or --buildings")
# One ray has no phase to be unsure of: both ends of the interval are the level.
add_test(NAME buildings_far_phase_lower_is_the_level COMMAND golden_compare integrated/buildings_far_phase_lower.csv integrated/buildings_far_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_far_phase_upper_is_the_level COMMAND golden_compare integrated/buildings_far_phase_upper.csv integrated/buildings_far_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_far_phase_lower_is_the_level buildings_far_phase_upper_is_the_level PROPERTIES FIXTURES_REQUIRED buildings_far_runs)
# Implicit interval filenames must not replace another metric or trajectory.
foreach(end lower upper)
  add_test(NAME buildings_phase_${end}_refuses_csel_collision COMMAND ${URBAN_RUN} --buildings tests/data/urban_far_building.csv --phase_interval 90 -o integrated/phase_collision.csv --csel_file ${CMAKE_CURRENT_LIST_DIR}/integrated/phase_collision_phase_${end}.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  add_test(NAME buildings_phase_${end}_refuses_trajectory_collision COMMAND ${URBAN_RUN} -t tracks/S-76D/197421AC.csv --buildings tests/data/urban_far_building.csv --phase_interval 90 -o integrated/phase_collision.csv -o integrated/phase_collision_phase_${end}.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(buildings_phase_${end}_refuses_csel_collision buildings_phase_${end}_refuses_trajectory_collision PROPERTIES PASS_REGULAR_EXPRESSION "--phase_interval output .* would overwrite another output")
endforeach()
add_test(NAME buildings_phase_refuses_metric_collision COMMAND ${URBAN_RUN} --buildings tests/data/urban_far_building.csv --phase_interval 90 -o integrated/phase_collision.csv --lamax_file integrated/phase_collision_phase_upper.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_phase_refuses_companion_collision COMMAND ${URBAN_RUN} -t tracks/S-76D/197421AC.csv --buildings tests/data/urban_far_building.csv --phase_interval 90 -o integrated/phase_collision -o integrated/phase_collision.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_phase_refuses_metric_collision buildings_phase_refuses_companion_collision PROPERTIES PASS_REGULAR_EXPRESSION "--phase_interval output .* would overwrite another output")
add_test(NAME buildings_refuses_coverage COMMAND ${URBAN_RUN} --buildings tests/data/urban_far_building.csv --coverage_file integrated/never_coverage.csv -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_coverage PROPERTIES PASS_REGULAR_EXPRESSION "--buildings does not combine with --coverage_file")
# The measured coherence through a run: with one ray any sum is that ray's level, so the far
# building's run is the open one again, at 1e-9; and over the test district the interval's ends
# bracket every microphone's SEL, some strictly.
add_test(NAME buildings_far_measured_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_far_building.csv --floor_absorption 1 --ray_coherence measured --state_tolerance 0 --no-gpu -o integrated/buildings_far_measured.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_measured_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_buildings.csv --ray_coherence measured --phase_interval 90 --no-gpu -o integrated/buildings_measured.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_far_measured_run buildings_measured_run PROPERTIES FIXTURES_SETUP buildings_measured_runs)
add_test(NAME buildings_far_measured_is_open COMMAND golden_compare integrated/buildings_far_measured.csv integrated/buildings_far_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_far_measured_is_open PROPERTIES FIXTURES_REQUIRED "buildings_measured_runs;buildings_far_runs")
add_test(NAME buildings_measured_interval_brackets COMMAND ${CMAKE_COMMAND} -DLEVEL=integrated/buildings_measured.csv -DLOWER=integrated/buildings_measured_phase_lower.csv -DUPPER=integrated/buildings_measured_phase_upper.csv -P tests/check_interval.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_measured_interval_brackets PROPERTIES FIXTURES_REQUIRED buildings_measured_runs)
add_test(NAME buildings_refuses_inside COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --frame_bearing 210 -b tests/data/urban_inside.csv -L warn --buildings tests/data/urban_buildings.csv -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_inside PROPERTIES PASS_REGULAR_EXPRESSION "observer 2's microphone stands inside building 8")
add_test(NAME buildings_refuses_urban_table COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --urban_table tables/urban_open_test.nc -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_urban_table PROPERTIES PASS_REGULAR_EXPRESSION "--buildings does not combine with --urban_table")
# --edge_source in the run: with one building 1000 km off nothing diffracts, so the edge source changes nothing
# (to the bit); it goes with --buildings, and not with --gpu.
add_test(NAME buildings_far_edge_source_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_far_building.csv --floor_absorption 1 --state_tolerance 0 --edge_source -o integrated/buildings_far_edge_source.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_far_edge_source_run PROPERTIES FIXTURES_SETUP buildings_far_edge_source_runs)
add_test(NAME buildings_far_edge_source_is_utd COMMAND golden_compare integrated/buildings_far_edge_source.csv integrated/buildings_far_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_far_edge_source_is_utd PROPERTIES FIXTURES_REQUIRED "buildings_far_edge_source_runs;buildings_far_runs")
# And with a building between the track and two microphones in its shadow it is wired through: the
# edge source moves their levels from UTD's (by ~1e-3 dB; golden_compare's mismatch is the pass).
add_test(NAME buildings_near_edge_source_run COMMAND ${URBAN_SHORT_BARE} --buildings tests/data/urban_near_building.csv --band_limits 500 2000 --edge_source -o integrated/buildings_near_edge_source.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_near_utd_run COMMAND ${URBAN_SHORT_BARE} --buildings tests/data/urban_near_building.csv --band_limits 500 2000 --no-gpu -o integrated/buildings_near_utd.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_near_edge_source_run buildings_near_utd_run PROPERTIES FIXTURES_SETUP buildings_near_edge_source_runs)
add_test(NAME buildings_edge_source_is_wired COMMAND golden_compare integrated/buildings_near_edge_source.csv integrated/buildings_near_utd.csv 1e-7 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_edge_source_is_wired PROPERTIES FIXTURES_REQUIRED buildings_near_edge_source_runs PASS_REGULAR_EXPRESSION "mismatch at value")
add_test(NAME buildings_refuses_lone_edge_source COMMAND ${URBAN_RUN} --edge_source -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_lone_edge_source PROPERTIES PASS_REGULAR_EXPRESSION "--edge_source, --state_tolerance, --path_radius, --gpu and --no-gpu go with --buildings")
add_test(NAME buildings_refuses_edge_source_gpu COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --edge_source --gpu -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_edge_source_gpu PROPERTIES PASS_REGULAR_EXPRESSION "--edge_source is the CPU's; --gpu does not take it")
# A ground model among the buildings (--ground_effect, point_paths.hpp FloorModel) replaces the floor's
# energy coefficient, and the rays' phases are the source's, not the ground model's turbulence.
add_test(NAME buildings_refuses_floor_absorption_with_ground COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --ground_effect --floor_absorption 0.2 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_floor_absorption_with_ground PROPERTIES PASS_REGULAR_EXPRESSION "--floor_absorption is the floor's energy coefficient among the buildings, which the ground model")
add_test(NAME buildings_refuses_ground_turbulence COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --ground_effect --turbulence 1e-4 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_ground_turbulence PROPERTIES PASS_REGULAR_EXPRESSION "does not take --turbulence or --gamma_t")
# With one building 1000 km off and the ground model's floor, each microphone hears the direct ray and the
# floor's, its Q at the sub-bands' middles, every pair's phases known: the open-ground run (--exact_ground)
# but for the floor's ray reading the sphere along its own launch and Q's change across a band, within
# 0.1 dB of SEL at the 4 ft microphones (0.060 measured), and the flush ones' plates within 0.003 dB.
# Over two ground classes, the 4 ft microphones on concrete too, within 0.15 dB (0.141).
set(URBAN_CLASS_RUN niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tests/data/urban_track_short.csv --frame_bearing 210 -b tests/data/urban_observers_classes.csv --ground_classes tests/data/ground_classes_two.csv ${URBAN_PLATE} -L warn)
add_test(NAME buildings_far_ground_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_far_building.csv --ground_effect --ray_coherence coherent --state_tolerance 0 --no-gpu -o integrated/buildings_far_ground.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_far_ground_reference COMMAND ${URBAN_DISTRICT_SHORT} --ground_effect --exact_ground -o integrated/buildings_far_ground_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_far_classes_run COMMAND ${URBAN_CLASS_RUN} --buildings tests/data/urban_far_building.csv --ray_coherence coherent --state_tolerance 0 --no-gpu -o integrated/buildings_far_classes.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_far_classes_reference COMMAND ${URBAN_CLASS_RUN} --exact_ground -o integrated/buildings_far_classes_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_far_ground_run buildings_far_ground_reference buildings_far_classes_run buildings_far_classes_reference PROPERTIES FIXTURES_SETUP buildings_far_ground_runs)
add_test(NAME buildings_far_ground_is_open COMMAND golden_compare integrated/buildings_far_ground.csv integrated/buildings_far_ground_reference.csv 0.1 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME buildings_far_classes_are_open COMMAND golden_compare integrated/buildings_far_classes.csv integrated/buildings_far_classes_reference.csv 0.15 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_far_ground_is_open buildings_far_classes_are_open PROPERTIES FIXTURES_REQUIRED buildings_far_ground_runs)
add_test(NAME buildings_refuses_lone_absorption COMMAND ${URBAN_RUN} --floor_absorption 0.2 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_lone_absorption PROPERTIES PASS_REGULAR_EXPRESSION "--wall_absorption and --floor_absorption go with --buildings")
add_test(NAME buildings_refuses_lone_coherence COMMAND ${URBAN_RUN} --ray_coherence energy -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_lone_coherence PROPERTIES PASS_REGULAR_EXPRESSION "--ray_coherence, --coherence_angle and --phase_interval go with --buildings")
add_test(NAME buildings_refuses_interval_with_known_phases COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --ray_coherence coherent --phase_interval 90 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_interval_with_known_phases PROPERTIES PASS_REGULAR_EXPRESSION "--phase_interval goes with --ray_coherence cone or measured")
add_test(NAME buildings_refuses_angle_without_cone COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --ray_coherence energy --coherence_angle 2 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_angle_without_cone PROPERTIES PASS_REGULAR_EXPRESSION "--coherence_angle goes with --ray_coherence cone")
add_test(NAME urban_buildings_refuses_inside COMMAND niceops_urban_table --buildings tests/data/urban_buildings.csv --receivers tests/data/urban_inside.csv --bands 1000 --frame_bearing 210 -o integrated/never.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings_refuses_inside PROPERTIES PASS_REGULAR_EXPRESSION "row 2: the microphone stands inside building 8")
# Receivers over and on roofs, flat and pitched (BuildingPathOptions::roofs): the table takes their
# reflections; without its roofs (--no_roofs) the builder refuses them.
add_test(NAME urban_buildings_roof_receivers_build COMMAND niceops_urban_table --buildings tests/data/urban_buildings.csv --roofs tests/data/urban_roofs.csv --receivers tests/data/urban_over_roofs.csv --bands 250 1000 --frame_bearing 210 -j 2 -o tables/urban_roof_receivers_test.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME urban_buildings_no_roofs_refuses_roof_receiver COMMAND niceops_urban_table --buildings tests/data/urban_buildings.csv --receivers tests/data/urban_over_roofs.csv --bands 1000 --frame_bearing 210 --no_roofs -o integrated/never.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings_no_roofs_refuses_roof_receiver PROPERTIES PASS_REGULAR_EXPRESSION "row 1: the microphone stands on or over the roof of building 1")
add_test(NAME urban_buildings_refuses_wave COMMAND niceops_urban_table --buildings tests/data/urban_buildings.csv --receivers tests/data/urban_observers.csv --bands 50 63 --wave_below 63 -o integrated/never.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings_refuses_wave PROPERTIES PASS_REGULAR_EXPRESSION "--wave_below solves ideal canyons only, not --buildings")
add_test(NAME urban_buildings_needs_receivers COMMAND niceops_urban_table --buildings tests/data/urban_buildings.csv --bands 1000 -o integrated/never.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_buildings_needs_receivers PROPERTIES PASS_REGULAR_EXPRESSION "--receivers goes with --buildings, and --buildings needs it")
# Refusals: a frame other than the table's, a run none of whose microphones it holds, the options
# stage 1 does not combine with, and a builder column in the wrong unit.
add_test(NAME urban_table_refuses_another_frame COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -b tests/data/urban_observers.csv --urban_table tables/urban_open_test.nc -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_refuses_another_frame PROPERTIES FIXTURES_REQUIRED urban_table_files PASS_REGULAR_EXPRESSION "was built with local \\+x along 210 degrees, and this run's \\+x is along 90 degrees \\(east, as no --frame_bearing was given\\)")
add_test(NAME urban_table_refuses_no_receivers COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --frame_bearing 210 --gridstep 4000 --urban_table tables/urban_open_test.nc -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_refuses_no_receivers PROPERTIES FIXTURES_REQUIRED urban_table_files PASS_REGULAR_EXPRESSION "none of its 7 receivers stands within 1 mm")
add_test(NAME urban_table_refuses_refraction COMMAND ${URBAN_RUN} --urban_table tables/urban_open_test.nc --refraction linear -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_refuses_refraction PROPERTIES PASS_REGULAR_EXPRESSION "--urban_table does not combine with --refraction yet")
add_test(NAME urban_table_refuses_db_cutoff COMMAND ${URBAN_RUN} --urban_table tables/urban_open_test.nc --db_cutoff 60 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_refuses_db_cutoff PROPERTIES PASS_REGULAR_EXPRESSION "--urban_table does not combine with --db_cutoff")
add_test(NAME urban_table_builder_checks_units COMMAND niceops_urban_table --canyons tests/data/urban_canyons.csv --bands 1000 --units m -o integrated/never.nc WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(urban_table_builder_checks_units PROPERTIES PASS_REGULAR_EXPRESSION "column 'x_ft' is in ft, and --units is m")
# niceops_point_field (tools/point_field.cpp): a point source in the test district's street, at
# its microphones on two threads with every ray written out, through the urban tools' shared
# district reader; and its refusal of a source inside a building.
add_test(NAME point_field_runs COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --source -1000 150 10 --microphones tests/data/urban_observers.csv --bands 125 250 500 1000 --rays integrated/point_field_rays.csv -j 2 --humidity 70 --no-gpu -o integrated/point_field.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME point_field_refuses_inside COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --source -1060 230 10 --microphones tests/data/urban_observers.csv --bands 125 --no_air_absorption -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_refuses_inside PROPERTIES PASS_REGULAR_EXPRESSION "--source stands inside building 4")
# The air's absorption depends strongly on the humidity, so a run without one is refused rather than assumed.
add_test(NAME point_field_refuses_no_humidity COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --source -1000 150 10 --microphones tests/data/urban_observers.csv --bands 125 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_refuses_no_humidity PROPERTIES PASS_REGULAR_EXPRESSION "needs --humidity")
# The edge-source integral runs on the CPU only, and is refused with --gpu
# whether or not the build has a GPU.
add_test(NAME point_field_edge_source_gpu_refused COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --source -1000 150 10 --microphones tests/data/urban_observers.csv --bands 125 --edge_source --gpu --no_air_absorption -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_edge_source_gpu_refused PROPERTIES PASS_REGULAR_EXPRESSION "does not run with --gpu or --fade")
# A source that moves (--sources): each row of tests/data/point_sources.csv over every microphone, the
# first the point_field_runs source.  The library test `a field from ends traced once is the pair's own`
# holds the fields to the single-source ones ray for ray.
# Scattering walls: the diffuse field lit by the source's specular paths (diffuse_facades.hpp).
add_test(NAME point_field_scattering_runs COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --source -1000 150 10 --microphones tests/data/urban_observers.csv --bands 125 1000 --scattering 0.3 --patch_size 20 -j 2 --humidity 70 --no-gpu -o integrated/point_field_scattering.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME point_field_scattering_out_of_range COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --source -1000 150 10 --microphones tests/data/urban_observers.csv --bands 125 --scattering 1.5 --no_air_absorption -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_scattering_out_of_range PROPERTIES PASS_REGULAR_EXPRESSION "--scattering must lie within")
add_test(NAME point_field_scattering_spectrum_runs COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --source -1000 150 10 --microphones tests/data/urban_observers.csv --bands 125 1000 --scattering_spectrum tests/data/scattering_spectrum.csv --patch_size 20 -j 2 --humidity 70 --no-gpu -o integrated/point_field_scattering_spectrum.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME point_field_scattering_spectrum_short COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --source -1000 150 10 --microphones tests/data/urban_observers.csv --bands 125 4000 --scattering_spectrum tests/data/scattering_spectrum.csv --no_air_absorption -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_scattering_spectrum_short PROPERTIES PASS_REGULAR_EXPRESSION "covers 100 to 2000 Hz, not the band at 4000 Hz")
# A wall beyond --diffuse_radius returns nothing diffusely, so it keeps its whole rho: with none in
# reach the levels are the specular run's.
add_test(NAME point_field_specular_reference COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --source -1000 150 10 --microphones tests/data/urban_observers.csv --bands 125 1000 -j 2 --humidity 70 --no-gpu -o integrated/point_field_specular_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME point_field_scattering_out_of_reach COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --source -1000 150 10 --microphones tests/data/urban_observers.csv --bands 125 1000 --scattering 0.5 --diffuse_radius 0.001 -j 2 --humidity 70 --no-gpu -o integrated/point_field_scattering_out_of_reach.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_specular_reference point_field_scattering_out_of_reach PROPERTIES FIXTURES_SETUP point_field_reach_runs)
add_test(NAME point_field_scattering_out_of_reach_is_specular COMMAND golden_compare integrated/point_field_scattering_out_of_reach.csv integrated/point_field_specular_reference.csv 1e-9 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_scattering_out_of_reach_is_specular PROPERTIES FIXTURES_REQUIRED point_field_reach_runs)
# niceops_fdtd_field (tools/fdtd_field.cpp): a block on the floor, a coarse grid on the CPU, half
# the space by the source's plane of symmetry.
add_test(NAME fdtd_field_runs COMMAND niceops_fdtd_field --buildings tests/data/fdtd_block.csv --microphones tests/data/fdtd_microphones.csv --source 0 0 0.1 --bands 250 500 --units m --cell 0.02 --margin 0.2 --mirror -j 2 --no-gpu -o integrated/fdtd_field.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME fdtd_field_order2_runs COMMAND niceops_fdtd_field --buildings tests/data/fdtd_block.csv --microphones tests/data/fdtd_microphones.csv --source 0 0 0.1 --bands 250 500 --units m --cell 0.02 --margin 0.2 --mirror --order 2 -j 2 --no-gpu -o integrated/fdtd_field_order2.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME fdtd_field_order_refused COMMAND niceops_fdtd_field --buildings tests/data/fdtd_block.csv --microphones tests/data/fdtd_microphones.csv --source 0 0 0.1 --bands 250 --units m --order 3 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(fdtd_field_order_refused PROPERTIES WILL_FAIL TRUE)
add_test(NAME point_field_sources_runs COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --sources tests/data/point_sources.csv --microphones tests/data/urban_observers.csv --bands 125 250 500 1000 --rays integrated/point_field_sources_rays.csv -j 2 --humidity 70 --no-gpu -o integrated/point_field_sources.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# The shared faces sampled (the default) against their integral at every wavenumber
# (--exact_shared_faces): within 0.1 dB in every band.
add_test(NAME point_field_sources_exact_runs COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --sources tests/data/point_sources.csv --microphones tests/data/urban_observers.csv --bands 125 250 500 1000 --exact_shared_faces -j 2 --humidity 70 --no-gpu -o integrated/point_field_sources_exact.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_sources_runs point_field_sources_exact_runs PROPERTIES FIXTURES_SETUP point_field_sources)
add_test(NAME point_field_shared_faces_sampled COMMAND golden_compare integrated/point_field_sources.csv integrated/point_field_sources_exact.csv 0.1 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_shared_faces_sampled PROPERTIES FIXTURES_REQUIRED point_field_sources)
# A relief that keeps everything in its specular order (tests/data/relief_flat.csv) changes nothing, in the
# run and in niceops_point_field; a relief's grating rays and specular shares on the GPU are the CPU's, within
# 0.01 dB of SEL; --facade_relief goes with --buildings and not with a scattering spectrum.
add_test(NAME point_field_relief_flat_runs COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --sources tests/data/point_sources.csv --microphones tests/data/urban_observers.csv --bands 125 250 500 1000 -j 2 --humidity 70 --no-gpu --relief tests/data/relief_flat.csv -o integrated/point_field_relief_flat.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_relief_flat_runs PROPERTIES FIXTURES_SETUP point_field_relief_flat)
add_test(NAME point_field_relief_flat_is_plain COMMAND golden_compare integrated/point_field_relief_flat.csv integrated/point_field_sources.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_relief_flat_is_plain PROPERTIES FIXTURES_REQUIRED "point_field_relief_flat;point_field_sources")
add_test(NAME buildings_relief_flat_run COMMAND ${URBAN_DISTRICT_SHORT} --buildings tests/data/urban_buildings.csv --facade_relief tests/data/relief_flat.csv --no-gpu -o integrated/buildings_relief_flat.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_relief_flat_run PROPERTIES FIXTURES_SETUP buildings_relief_flat_runs)
add_test(NAME buildings_relief_flat_is_plain COMMAND golden_compare integrated/buildings_relief_flat.csv integrated/buildings_short.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_relief_flat_is_plain PROPERTIES FIXTURES_REQUIRED "buildings_relief_flat_runs;buildings_short_runs")
add_test(NAME buildings_refuses_lone_relief COMMAND ${URBAN_RUN} --facade_relief tests/data/relief_test.csv -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_lone_relief PROPERTIES PASS_REGULAR_EXPRESSION "--facade_relief, --relief_scattering, --grating_reflections, --patch_size and --diffuse_radius go with --buildings")
add_test(NAME buildings_refuses_relief_and_spectrum COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --facade_relief tests/data/relief_test.csv --facade_scattering_spectrum tests/data/scattering_spectrum.csv -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(buildings_refuses_relief_and_spectrum PROPERTIES PASS_REGULAR_EXPRESSION "--facade_relief and --facade_scattering_spectrum are each the walls' scattering")
# A floor of finite impedance (--flow_resistance) whose flow resistance is enormous reflects as the rigid
# floor does, to 0.01 dB in every band (0.006 measured), its bounces' Q all but 1; with
# --floor_absorption it is refused.
add_test(NAME point_field_hard_floor_runs COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --sources tests/data/point_sources.csv --microphones tests/data/urban_observers.csv --bands 125 250 500 1000 -j 2 --humidity 70 --no-gpu --flow_resistance 1e8 -o integrated/point_field_hard_floor.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_hard_floor_runs PROPERTIES FIXTURES_SETUP point_field_hard_floor)
add_test(NAME point_field_hard_floor_is_rigid COMMAND golden_compare integrated/point_field_hard_floor.csv integrated/point_field_sources.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_hard_floor_is_rigid PROPERTIES FIXTURES_REQUIRED "point_field_hard_floor;point_field_sources")
add_test(NAME point_field_refuses_two_floors COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --source 0 0 10 --microphones tests/data/urban_observers.csv --bands 125 --no_air_absorption --flow_resistance 200 --floor_absorption 0.1 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_refuses_two_floors PROPERTIES PASS_REGULAR_EXPRESSION "--flow_resistance gives the floor an impedance")
add_test(NAME point_field_refuses_two_sources COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --source -1000 150 10 --sources tests/data/point_sources.csv --microphones tests/data/urban_observers.csv --bands 125 --no_air_absorption -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_refuses_two_sources PROPERTIES PASS_REGULAR_EXPRESSION "excludes")
add_test(NAME point_field_refuses_nonfinite COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --sources tests/data/point_sources_nan.csv --microphones tests/data/urban_observers.csv --bands 125 --no_air_absorption -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(point_field_refuses_nonfinite PROPERTIES PASS_REGULAR_EXPRESSION "row 2's y is not a finite number")
# --gpu: the same run on the GPU writes the same rays, and levels within
# 0.05 dB of the CPU's, its shared faces sampled as the CPU samples them
# (the GPU integrates the low bands over every pair,
# where the CPU's pair sum drops the pairs its sinc damps below 1/30; FIDELITY
# has the difference).  A build without a GPU backend refuses the flag.
if(NICEOPS_GPU_BACKEND STREQUAL "OFF")
  add_test(NAME point_field_gpu_refused COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --sources tests/data/point_sources.csv --microphones tests/data/urban_observers.csv --bands 125 --gpu --no_air_absorption -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(point_field_gpu_refused PROPERTIES PASS_REGULAR_EXPRESSION "this build has no GPU backend")
  add_test(NAME buildings_gpu_refused COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --gpu -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(buildings_gpu_refused PROPERTIES PASS_REGULAR_EXPRESSION "this build has no GPU backend")
else()
  add_test(NAME point_field_gpu_runs COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --sources tests/data/point_sources.csv --microphones tests/data/urban_observers.csv --bands 125 250 500 1000 --rays integrated/point_field_sources_gpu_rays.csv -j 2 --gpu --humidity 70 -o integrated/point_field_sources_gpu.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(point_field_gpu_runs PROPERTIES FIXTURES_SETUP point_field_gpu)
  add_test(NAME point_field_gpu_matches_cpu COMMAND golden_compare integrated/point_field_sources_gpu.csv integrated/point_field_sources.csv 0.05 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  add_test(NAME point_field_gpu_same_rays COMMAND ${CMAKE_COMMAND} -E compare_files integrated/point_field_sources_gpu_rays.csv integrated/point_field_sources_rays.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(point_field_gpu_matches_cpu point_field_gpu_same_rays PROPERTIES FIXTURES_REQUIRED "point_field_gpu;point_field_sources")
  # With no diffraction a batch has no wedges, whose buffer every dispatch of
  # the amplitudes still binds.
  add_test(NAME point_field_gpu_no_diffraction_runs COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --sources tests/data/point_sources.csv --microphones tests/data/urban_observers.csv --bands 125 250 --no_diffraction -j 2 --gpu --humidity 70 -o integrated/point_field_specular_gpu.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  add_test(NAME point_field_no_diffraction_runs COMMAND niceops_point_field --buildings tests/data/urban_buildings.csv --sources tests/data/point_sources.csv --microphones tests/data/urban_observers.csv --bands 125 250 --no_diffraction -j 2 --humidity 70 --no-gpu -o integrated/point_field_specular.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(point_field_gpu_no_diffraction_runs point_field_no_diffraction_runs PROPERTIES FIXTURES_SETUP point_field_specular)
  add_test(NAME point_field_gpu_no_diffraction_matches_cpu COMMAND golden_compare integrated/point_field_specular_gpu.csv integrated/point_field_specular.csv 0.05 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(point_field_gpu_no_diffraction_matches_cpu PROPERTIES FIXTURES_REQUIRED point_field_specular)
  # niceops --buildings on the GPU (the default where there is one) against the CPU (--no-gpu):
  # the demo approach's SELs within 0.01 dB; and with every state evaluated, within 0.01 dB too.
  add_test(NAME buildings_gpu_run COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --gpu -o integrated/buildings_gpu.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  add_test(NAME buildings_gpu_every_state_run COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --gpu --state_tolerance 0 -o integrated/buildings_gpu_every_state.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(buildings_gpu_run buildings_gpu_every_state_run PROPERTIES FIXTURES_SETUP buildings_gpu_runs)
  add_test(NAME buildings_gpu_matches_cpu COMMAND golden_compare integrated/buildings_gpu.csv integrated/buildings_sampled.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  add_test(NAME buildings_gpu_every_state_matches_cpu COMMAND golden_compare integrated/buildings_gpu_every_state.csv integrated/buildings_every_state.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(buildings_gpu_matches_cpu buildings_gpu_every_state_matches_cpu PROPERTIES FIXTURES_REQUIRED "buildings_gpu_runs;buildings_sampling_runs")
  # With a relief scattering by the hybrid (grating rays and the diffuse field): the GPU's SELs are the CPU's.
  add_test(NAME buildings_relief_gpu_run COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --facade_relief tests/data/relief_test.csv --gpu -o integrated/buildings_relief_gpu.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  add_test(NAME buildings_relief_cpu_run COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --facade_relief tests/data/relief_test.csv --no-gpu -o integrated/buildings_relief_cpu.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(buildings_relief_gpu_run buildings_relief_cpu_run PROPERTIES FIXTURES_SETUP buildings_relief_gpu_runs)
  add_test(NAME buildings_relief_gpu_matches_cpu COMMAND golden_compare integrated/buildings_relief_gpu.csv integrated/buildings_relief_cpu.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(buildings_relief_gpu_matches_cpu PROPERTIES FIXTURES_REQUIRED buildings_relief_gpu_runs)
  # And with the ground model's floor over two classes, the flush microphones on plates: the GPU takes
  # each ray's Q and P_d or P_r per sub-band from the CPU, within 0.01 dB of SEL (1e-4 measured).
  add_test(NAME buildings_ground_gpu_run COMMAND ${URBAN_CLASS_RUN} --buildings tests/data/urban_buildings.csv --gpu -o integrated/buildings_ground_gpu.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  add_test(NAME buildings_ground_cpu_run COMMAND ${URBAN_CLASS_RUN} --buildings tests/data/urban_buildings.csv --no-gpu -o integrated/buildings_ground_cpu.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(buildings_ground_gpu_run buildings_ground_cpu_run PROPERTIES FIXTURES_SETUP buildings_ground_gpu_runs)
  add_test(NAME buildings_ground_gpu_matches_cpu COMMAND golden_compare integrated/buildings_ground_gpu.csv integrated/buildings_ground_cpu.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(buildings_ground_gpu_matches_cpu PROPERTIES FIXTURES_REQUIRED buildings_ground_gpu_runs)
  # The measured coherence and the phase interval are the CPU's: --gpu asked for is refused there,
  # as the arguments are read rather than track by track.
  add_test(NAME buildings_gpu_refuses_measured COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --ray_coherence measured --gpu -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  add_test(NAME buildings_gpu_refuses_interval COMMAND ${URBAN_DISTRICT} --buildings tests/data/urban_buildings.csv --phase_interval 90 --gpu -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(buildings_gpu_refuses_measured buildings_gpu_refuses_interval PROPERTIES PASS_REGULAR_EXPRESSION "--gpu: the GPU does not take --ray_coherence measured or --phase_interval" FAIL_REGULAR_EXPRESSION "Skipping")
endif()
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
# A WGS84 track with a heading column: the geodetic demo track with its own ground-track bearing
# (central differences of its positions) and that turned 15 degrees.  Along the track, heading
# mode is the default nose to within the finite-difference bearing (median 0.04 dB, 0.55 at worst
# on this grid); crabbed, it moves the footprint by up to 2.9 dB.
set(GEODETIC_FRAME --origin 39.82914323 -74.40996141 --frame_bearing 210 --ground_elevation 54.2651 --gridstep 328 -L warn)
add_test(NAME geodetic_heading_default COMMAND niceops -t tests/data/demo_geodetic_heading.csv ${GEODETIC_FRAME} -o integrated/geodetic_heading_default.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME geodetic_heading_on_track COMMAND niceops -t tests/data/demo_geodetic_heading.csv ${GEODETIC_FRAME} --nose heading -o integrated/geodetic_heading_on_track.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME geodetic_heading_crabbed COMMAND niceops -t tests/data/demo_geodetic_heading_crabbed.csv ${GEODETIC_FRAME} --nose heading -o integrated/geodetic_heading_crabbed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(geodetic_heading_default geodetic_heading_on_track geodetic_heading_crabbed PROPERTIES FIXTURES_SETUP geodetic_heading_runs)
add_test(NAME geodetic_heading_on_track_is_velocity COMMAND golden_compare integrated/geodetic_heading_on_track.csv integrated/geodetic_heading_default.csv 0.6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME geodetic_heading_crabbed_differs COMMAND golden_compare integrated/geodetic_heading_crabbed.csv integrated/geodetic_heading_default.csv 0.6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(geodetic_heading_crabbed_differs PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value")
set_tests_properties(geodetic_heading_on_track_is_velocity geodetic_heading_crabbed_differs PROPERTIES FIXTURES_REQUIRED geodetic_heading_runs)
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
# naming it, and the default one-sided law, must reproduce the demo run (on the coarse grid); nearest and aam must run and
# move it (golden_compare's mismatch is the pass), and differ from each other.  The demo database is
# S-76D_M3, whose thrust axis is spheres, so aam there searches its thrust layers.
add_test(NAME blend_barycentric COMMAND niceops --no-spectral --no-ground_effect --blend barycentric --load_factor_correction one_sided ${DEMO_COARSE} -L warn -o integrated/blend_barycentric.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME blend_nearest COMMAND niceops --no-spectral --no-ground_effect --blend nearest ${DEMO_COARSE} -L warn -o integrated/blend_nearest.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME blend_aam COMMAND niceops --no-spectral --no-ground_effect --blend aam ${DEMO_COARSE} -L warn -o integrated/blend_aam.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(blend_barycentric blend_nearest blend_aam PROPERTIES FIXTURES_SETUP blend_runs)
add_test(NAME blend_barycentric_is_golden COMMAND golden_compare integrated/blend_barycentric.csv integrated/demo_coarse.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME blend_nearest_moves_footprint COMMAND golden_compare integrated/blend_nearest.csv integrated/demo_coarse.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME blend_aam_moves_footprint COMMAND golden_compare integrated/blend_aam.csv integrated/demo_coarse.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME blend_aam_is_not_nearest COMMAND golden_compare integrated/blend_aam.csv integrated/blend_nearest.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(blend_nearest_moves_footprint blend_aam_moves_footprint blend_aam_is_not_nearest PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value")
set_tests_properties(blend_barycentric_is_golden blend_nearest_moves_footprint blend_aam_moves_footprint blend_aam_is_not_nearest PROPERTIES FIXTURES_REQUIRED "blend_runs;demo_coarse_run")
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
# --ground_classes (src/ground_classes.hpp) and the observer file's sixth column
# (id,x,y,z,height,ground) name each other; the reader's own refusals are unit
# tests ([ground_classes]).
set(GROUND_CLASS_DEMO -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv -L warn)
set(GROUND_CLASS_GRASS --impedance_model variable_porosity --flow_resistance 200)
add_test(NAME ground_classes_column_needs_file COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} -b tests/data/ground_classes_mixed.csv -o integrated/ground_classes_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_column_needs_file PROPERTIES PASS_REGULAR_EXPRESSION "gives each microphone a ground class .* and there is no --ground_classes")
add_test(NAME ground_classes_refuses_grid COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --ground_classes tests/data/ground_classes_two.csv --gridstep 4000 -o integrated/ground_classes_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_refuses_grid PROPERTIES PASS_REGULAR_EXPRESSION "grid observers have none")
add_test(NAME ground_classes_need_column COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --ground_classes tests/data/ground_classes_two.csv -b tests/data/mixed_receivers.csv -o integrated/ground_classes_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_need_column PROPERTIES PASS_REGULAR_EXPRESSION "--ground_classes needs the observer file to name each microphone's class in a sixth column")
add_test(NAME ground_classes_refuse_undefined_class COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --ground_classes tests/data/ground_classes_one.csv -b tests/data/ground_classes_mixed.csv -o integrated/ground_classes_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_refuse_undefined_class PROPERTIES PASS_REGULAR_EXPRESSION "data row 1 stands on ground class 0, which --ground_classes tests/data/ground_classes_one.csv does not define \\(it defines 7\\)")
add_test(NAME ground_classes_refuse_fractional_class COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --ground_classes tests/data/ground_classes_two.csv -b tests/data/ground_classes_bad_id.csv -o integrated/ground_classes_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_refuse_fractional_class PROPERTIES PASS_REGULAR_EXPRESSION "data row 2 gives the ground class 2.5; a class is an integer")
# Each microphone over its class.  One class repeating the command line's ground is the run
# without classes, byte for byte: inherited (ground_classes_one.csv) and with every column given
# at the command line's values in meters (ground_classes_one_explicit.csv; under --units ft the
# two lengths go through x * 0.3048 and need not be), for point receivers with their tables,
# for plates, and over the hills.
add_test(NAME ground_class_scalar_point COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} -b tests/data/mixed_receivers.csv -o integrated/ground_class_scalar_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ground_class_single_point COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --ground_classes tests/data/ground_classes_one.csv -b tests/data/mixed_receivers_class7.csv -o integrated/ground_class_single_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ground_class_explicit_point COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --ground_classes tests/data/ground_classes_one_explicit.csv -b tests/data/mixed_receivers_class7.csv -o integrated/ground_class_explicit_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ground_class_scalar_plate COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --plate_table plates/gr1425_vp200.nc -b tests/data/mixed_receivers.csv -o integrated/ground_class_scalar_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ground_class_single_plate COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --plate_table plates/gr1425_vp200.nc --ground_classes tests/data/ground_classes_one.csv -b tests/data/mixed_receivers_class7.csv -o integrated/ground_class_single_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_class_scalar_point ground_class_single_point ground_class_explicit_point ground_class_scalar_plate ground_class_single_plate PROPERTIES FIXTURES_SETUP ground_class_single_runs)
add_test(NAME ground_class_single_is_scalar_point COMMAND ${CMAKE_COMMAND} -E compare_files integrated/ground_class_single_point.csv integrated/ground_class_scalar_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ground_class_explicit_is_scalar_point COMMAND ${CMAKE_COMMAND} -E compare_files integrated/ground_class_explicit_point.csv integrated/ground_class_scalar_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ground_class_single_is_scalar_plate COMMAND ${CMAKE_COMMAND} -E compare_files integrated/ground_class_single_plate.csv integrated/ground_class_scalar_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_class_single_is_scalar_point ground_class_explicit_is_scalar_point ground_class_single_is_scalar_plate PROPERTIES FIXTURES_REQUIRED ground_class_single_runs)
# A grass and concrete array (ground_classes_two.csv: class 0 inherits the command line's grass,
# class 1 is Delany-Bazley at 10000): plates and 4 ft poles on grass, then 4 ft and 10 ft poles
# on concrete, so one height stands on both classes.  With --exact_ground it has to equal a grass
# run over the grass rows and a concrete run over the concrete rows; with its tables it has to
# agree with that to the tables' accuracy, which a table looked up by height alone would miss by
# decibels at the concrete 4 ft poles.  It must differ from the all-grass run first at the first
# concrete microphone's level (value 35), which an admittance taken from the wrong class would
# not.  And its peak search, bounded over both classes, must equal the march, flat and over hills.
add_test(NAME ground_classes_mixed COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --plate_table plates/gr1425_vp200.nc --ground_classes tests/data/ground_classes_two.csv -b tests/data/ground_classes_mixed.csv -o integrated/ground_classes_mixed.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ground_classes_mixed_marched COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --plate_table plates/gr1425_vp200.nc --ground_classes tests/data/ground_classes_two.csv -b tests/data/ground_classes_mixed.csv -D integrated/ground_classes_mixed_dnl.csv -o integrated/ground_classes_mixed_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ground_classes_mixed_exact COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --exact_ground --plate_table plates/gr1425_vp200.nc --ground_classes tests/data/ground_classes_two.csv -b tests/data/ground_classes_mixed.csv -o integrated/ground_classes_mixed_exact.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ground_classes_grass_only COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --exact_ground --plate_table plates/gr1425_vp200.nc -b tests/data/ground_classes_grass.csv -o integrated/ground_classes_grass.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ground_classes_concrete_only COMMAND niceops ${GROUND_CLASS_DEMO} --flow_resistance 10000 --exact_ground -b tests/data/ground_classes_concrete.csv -o integrated/ground_classes_concrete.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_mixed ground_classes_mixed_marched ground_classes_mixed_exact ground_classes_grass_only ground_classes_concrete_only PROPERTIES FIXTURES_SETUP ground_classes_runs)
add_test(NAME ground_classes_expected COMMAND ${CMAKE_COMMAND} -DOUT=integrated/ground_classes_expected.csv "-DINPUTS=integrated/ground_classes_grass.csv;integrated/ground_classes_concrete.csv" -P tests/concatenate.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_expected PROPERTIES FIXTURES_REQUIRED ground_classes_runs FIXTURES_SETUP ground_classes_expected_file)
add_test(NAME ground_classes_mixed_are_their_subsets COMMAND golden_compare integrated/ground_classes_mixed_exact.csv integrated/ground_classes_expected.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_mixed_are_their_subsets PROPERTIES FIXTURES_REQUIRED "ground_classes_runs;ground_classes_expected_file")
add_test(NAME ground_classes_tabled_matches_exact COMMAND golden_compare integrated/ground_classes_mixed.csv integrated/ground_classes_mixed_exact.csv 0.02 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ground_classes_windowed_is_exact COMMAND golden_compare integrated/ground_classes_mixed.csv integrated/ground_classes_mixed_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_tabled_matches_exact ground_classes_windowed_is_exact PROPERTIES FIXTURES_REQUIRED ground_classes_runs)
add_test(NAME ground_classes_impedance_is_the_observers COMMAND golden_compare integrated/ground_classes_mixed_exact.csv integrated/mixed_receivers.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_impedance_is_the_observers PROPERTIES FIXTURES_REQUIRED "ground_classes_runs;mixed_runs" PASS_REGULAR_EXPRESSION "mismatch at value 35:")
# What the classes cannot serve, each refused with its reason.
add_test(NAME ground_classes_refuse_refraction COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --ground_classes tests/data/ground_classes_two.csv -b tests/data/ground_classes_mixed.csv --atmosphere atmospheres/B407_uniform.csv --refraction linear -o integrated/ground_classes_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_refuse_refraction PROPERTIES PASS_REGULAR_EXPRESSION "--ground_classes does not combine with --refraction")
add_test(NAME ground_classes_refuse_table COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --ground_classes tests/data/ground_classes_two.csv -b tests/data/ground_classes_mixed.csv --propagation_table tables/absent.nc -o integrated/ground_classes_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_refuse_table PROPERTIES PASS_REGULAR_EXPRESSION "--ground_classes does not combine with --propagation_table: a table was solved over one ground")
add_test(NAME ground_classes_refuse_broadband COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --ground_classes tests/data/ground_classes_two.csv -b tests/data/ground_classes_mixed.csv --no-spectral -o integrated/ground_classes_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_refuse_broadband PROPERTIES PASS_REGULAR_EXPRESSION "--ground_classes needs the spectral path")
add_test(NAME ground_classes_refuse_no_ground COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --ground_classes tests/data/ground_classes_two.csv -b tests/data/ground_classes_mixed.csv --no-ground_effect -o integrated/ground_classes_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_refuse_no_ground PROPERTIES PASS_REGULAR_EXPRESSION "--no-ground_effect turns it off")
add_test(NAME ground_classes_need_spectra COMMAND niceops --humidity 20 --ground_classes tests/data/ground_classes_two.csv -b tests/data/ground_classes_mixed.csv -L warn -o integrated/ground_classes_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_need_spectra PROPERTIES PASS_REGULAR_EXPRESSION "--ground_classes needs a database carrying per-band spectra")
add_test(NAME ground_classes_plates_on_one_class COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --plate_table plates/gr1425_vp200.nc --ground_classes tests/data/ground_classes_two.csv -b tests/data/ground_classes_two_plates.csv -o integrated/ground_classes_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_plates_on_one_class PROPERTIES PASS_REGULAR_EXPRESSION "the flush microphones stand on ground classes 0 \\(grass\\) and 1 \\(concrete\\)")
add_test(NAME ground_classes_plate_matches_class COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --plate_table plates/gr1425_db225.nc --ground_classes tests/data/ground_classes_two.csv -b tests/data/ground_classes_mixed.csv -o integrated/ground_classes_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_plate_matches_class PROPERTIES PASS_REGULAR_EXPRESSION "--plate_table over ground class 0 \\(grass\\), which the flush microphones stand on: Plate table plates/gr1425_db225.nc was computed over")
add_test(NAME ground_classes_refuse_not_passive COMMAND niceops ${GROUND_CLASS_DEMO} ${GROUND_CLASS_GRASS} --ground_classes tests/data/ground_classes_not_passive.csv -b tests/data/ground_classes_mixed.csv -o integrated/ground_classes_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_refuse_not_passive PROPERTIES PASS_REGULAR_EXPRESSION "Ground class 1 \\(thin, tests/data/ground_classes_not_passive.csv line 3\\): The ground impedance model is not passive")
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
set_tests_properties(refraction_uniform_point PROPERTIES FIXTURES_SETUP refraction_point_runs)
add_test(NAME refraction_uniform_point_is_straight COMMAND golden_compare integrated/refraction_uniform_point.csv integrated/atmosphere_homogeneous.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(refraction_uniform_point_is_straight PROPERTIES FIXTURES_REQUIRED "refraction_point_runs;atmosphere_runs")
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
add_test(NAME stratified_uniform_point_is_straight COMMAND golden_compare integrated/stratified_uniform_point.csv integrated/atmosphere_homogeneous.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_uniform_point_is_straight PROPERTIES FIXTURES_REQUIRED "stratified_runs;atmosphere_runs")
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
add_test(NAME ray_stride_needs_stratified COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction linear --ray_stride 4 --gridstep 4000 -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ray_stride_needs_stratified PROPERTIES PASS_REGULAR_EXPRESSION "--ray_stride needs --refraction stratified")
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
set_tests_properties(stratified_resolved_point stratified_resolved_point_marched stratified_resolved_plate stratified_resolved_plate_marched stratified_resolved_straight_point stratified_resolved_straight_point_marched PROPERTIES FIXTURES_SETUP stratified_resolved_runs)
add_test(NAME stratified_resolved_point_windowed_is_exact COMMAND golden_compare integrated/stratified_resolved_point.csv integrated/stratified_resolved_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_resolved_plate_windowed_is_exact COMMAND golden_compare integrated/stratified_resolved_plate.csv integrated/stratified_resolved_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_resolved_straight_point_windowed_is_exact COMMAND golden_compare integrated/stratified_resolved_straight_point.csv integrated/stratified_resolved_straight_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_resolved_point_windowed_is_exact stratified_resolved_plate_windowed_is_exact stratified_resolved_straight_point_windowed_is_exact PROPERTIES FIXTURES_REQUIRED stratified_resolved_runs)
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
set_tests_properties(stratified_windy_plate stratified_calm_floor_plate PROPERTIES FIXTURES_SETUP stratified_calm_floor_runs)
add_test(NAME stratified_calm_floor_plate_is_floorless COMMAND golden_compare integrated/stratified_calm_floor_plate.csv integrated/stratified_windy_plate.csv 0.3 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_calm_floor_plate_is_floorless PROPERTIES FIXTURES_REQUIRED stratified_calm_floor_runs)
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
# --ray_stride is 4 by default.  Every state solved (--ray_stride 1) stays pinned by the goldens made
# before it was (*_exact.csv); the default is within 0.005 dB of it on the windy profile and 0.02 in the
# duct (measured 0.001 and 0.005), and not equal to it.
add_test(NAME stratified_windy_point_exact COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy.csv --refraction stratified --frame_bearing 210 --ground_effect --gridstep 984 -L warn --ray_stride 1 -o integrated/stratified_windy_point_exact.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_duct_point_exact COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_inversion.csv --refraction stratified --ground_effect --gridstep 984 -L warn --ray_stride 1 -o integrated/stratified_duct_point_exact.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_windy_point_exact stratified_duct_point_exact PROPERTIES FIXTURES_SETUP stratified_exact_runs)
add_test(NAME stratified_windy_point_exact_golden COMMAND golden_compare integrated/stratified_windy_point_exact.csv golden/stratified_windy_point_exact.csv 1e-4 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_duct_point_exact_golden COMMAND golden_compare integrated/stratified_duct_point_exact.csv golden/stratified_duct_point_exact.csv 1e-4 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_stride_windy_is_close COMMAND golden_compare integrated/stratified_windy_point.csv integrated/stratified_windy_point_exact.csv 0.005 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_stride_duct_is_close COMMAND golden_compare integrated/stratified_duct_point.csv integrated/stratified_duct_point_exact.csv 0.02 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME stratified_stride_is_not_exact COMMAND golden_compare integrated/stratified_windy_point.csv integrated/stratified_windy_point_exact.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(stratified_stride_is_not_exact PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value")
set_tests_properties(stratified_windy_point_exact_golden stratified_duct_point_exact_golden stratified_stride_windy_is_close stratified_stride_duct_is_close stratified_stride_is_not_exact PROPERTIES FIXTURES_REQUIRED "stratified_exact_runs;stratified_runs;stratified_duct_runs")

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
# niceops_build_atmosphere (tools/build_atmosphere.cpp).  One surface observation with --wind's
# choices is the profile wind_profile_written writes, byte for byte.  The harness layout (seven
# LIDAR gates, one a fill value, a string every 10 ft, a station humidity, the harness's level list)
# is the validation harness's own Weather.profile, run by tests/data/atmosphere_builder_expected.py;
# the same build in meters is the same profile.  A built profile runs, and the peak search over it
# is exact.
set(BUILD_ATMOSPHERE_DATA tests/data/atmosphere_builder)
set(BUILD_ATMOSPHERE_HARNESS build --column ${BUILD_ATMOSPHERE_DATA}/lidar_ft.csv --column ${BUILD_ATMOSPHERE_DATA}/string_ft.csv --missing 9999 --humidity 62.5 --levels 0,1,3,10,20 --z0 0.1 --top 5000)
add_test(NAME build_atmosphere_surface COMMAND niceops_build_atmosphere --units ft build --wind 12kt,250,6.5 --z0 0.1 --surface_layer_depth 3000 --levels_per_decade 10 --temperature 288.15 --humidity 40 --pressure 100 --above lapse --humidity_aloft held_concentration --top 3000 -o integrated/build_atmosphere_surface.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME build_atmosphere_harness COMMAND niceops_build_atmosphere ${BUILD_ATMOSPHERE_HARNESS} -o integrated/build_atmosphere_harness.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME build_atmosphere_harness_m COMMAND niceops_build_atmosphere --units m build --column ${BUILD_ATMOSPHERE_DATA}/lidar_m.csv --column ${BUILD_ATMOSPHERE_DATA}/string_m.csv --missing 9999 --humidity 62.5 --levels 0,0.3048,0.9144000000000001,3.048,6.096 --z0 0.030480000000000004 --top 1524 -o integrated/build_atmosphere_harness_m.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_surface build_atmosphere_harness build_atmosphere_harness_m PROPERTIES FIXTURES_SETUP build_atmosphere_runs)
add_test(NAME build_atmosphere_surface_is_wind COMMAND ${CMAKE_COMMAND} -E compare_files integrated/build_atmosphere_surface.csv integrated/wind_profile_law.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_surface_is_wind PROPERTIES FIXTURES_REQUIRED "build_atmosphere_runs;wind_law_runs")
add_test(NAME build_atmosphere_is_harness_profile COMMAND golden_compare integrated/build_atmosphere_harness.csv tests/data/atmosphere_builder_expected_ft.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME build_atmosphere_units_agree COMMAND golden_compare integrated/build_atmosphere_harness_m.csv tests/data/atmosphere_builder_expected_m.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set(BUILD_ATMOSPHERE_RUN -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere integrated/build_atmosphere_harness.csv --refraction stratified --frame_bearing 210 --ground_effect --gridstep 984 -L warn)
add_test(NAME build_atmosphere_stratified COMMAND niceops ${BUILD_ATMOSPHERE_RUN} --manifest integrated/build_atmosphere_manifest.json -o integrated/build_atmosphere_stratified.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME build_atmosphere_stratified_marched COMMAND niceops ${BUILD_ATMOSPHERE_RUN} -D integrated/build_atmosphere_stratified_dnl.csv -o integrated/build_atmosphere_stratified_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_stratified build_atmosphere_stratified_marched PROPERTIES FIXTURES_SETUP build_atmosphere_stratified_runs FIXTURES_REQUIRED build_atmosphere_runs)
set_tests_properties(build_atmosphere_stratified PROPERTIES FAIL_REGULAR_EXPRESSION "is calm at;no level at or below the ground")
add_test(NAME build_atmosphere_stratified_windowed_is_exact COMMAND golden_compare integrated/build_atmosphere_stratified.csv integrated/build_atmosphere_stratified_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_is_harness_profile build_atmosphere_units_agree PROPERTIES FIXTURES_REQUIRED build_atmosphere_runs)
set_tests_properties(build_atmosphere_stratified_windowed_is_exact PROPERTIES FIXTURES_REQUIRED build_atmosphere_stratified_runs)
# --atmosphere warns, without refusing, on what it takes on trust: no level at the ground (the
# lowest is held down), and a calm level above the ground under wind.  A built profile has
# neither (build_atmosphere_stratified, above, carries the FAIL_REGULAR_EXPRESSION); its sidecar is
# recorded in that run's manifest.
add_test(NAME atmosphere_warns_no_ground_level COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere tests/data/atmosphere_no_ground_ft.csv --gridstep 1640 -L warn -o integrated/atmosphere_no_ground.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(atmosphere_warns_no_ground_level PROPERTIES PASS_REGULAR_EXPRESSION "has no level at or below the ground: its lowest, at 50 ft, is held down to it")
add_test(NAME atmosphere_warns_calm_floor COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_windy_calm_floor.csv --gridstep 1640 -L warn -o integrated/atmosphere_calm_floor.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(atmosphere_warns_calm_floor PROPERTIES PASS_REGULAR_EXPRESSION "example_windy_calm_floor.csv is calm at 0.1 ft, above the ground")
add_test(NAME atmosphere_still_does_not_warn COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --atmosphere atmospheres/example_inversion.csv --gridstep 1640 -L warn -o integrated/atmosphere_still.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(atmosphere_still_does_not_warn PROPERTIES FAIL_REGULAR_EXPRESSION "is calm at;no level at or below the ground")
add_test(NAME build_atmosphere_manifest_records_provenance COMMAND ${CMAKE_COMMAND} -E cat integrated/build_atmosphere_manifest.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_manifest_records_provenance PROPERTIES FIXTURES_REQUIRED build_atmosphere_stratified_runs PASS_REGULAR_EXPRESSION "\"atmosphere_provenance\"")
# fit_z0 gives back the roughness length of runs that follow the log law exactly (written by
# tests/data/atmosphere_builder/fit_runs.py), leaving out the one in light wind.
add_test(NAME fit_z0_known COMMAND niceops_build_atmosphere fit_z0 ${BUILD_ATMOSPHERE_DATA}/fit_run1.csv ${BUILD_ATMOSPHERE_DATA}/fit_run2.csv ${BUILD_ATMOSPHERE_DATA}/fit_run3.csv ${BUILD_ATMOSPHERE_DATA}/fit_run4.csv ${BUILD_ATMOSPHERE_DATA}/fit_run5.csv ${BUILD_ATMOSPHERE_DATA}/fit_run6.csv --missing 9999 --draws 200 --table integrated/fit_z0_known.csv -o integrated/fit_z0_known.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(fit_z0_known PROPERTIES FIXTURES_SETUP fit_z0_runs PASS_REGULAR_EXPRESSION "excluded fit_run6: light wind")
add_test(NAME fit_z0_recovers_known COMMAND golden_compare integrated/fit_z0_known.csv ${BUILD_ATMOSPHERE_DATA}/fit_expected.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(fit_z0_recovers_known PROPERTIES FIXTURES_REQUIRED fit_z0_runs)
# --z0_from takes the fit's z0_m in meters whatever --units says: the harness build in feet from
# that fit (0.1 ft to 2e-13) is the harness profile.
add_test(NAME build_atmosphere_z0_from COMMAND niceops_build_atmosphere build --column ${BUILD_ATMOSPHERE_DATA}/lidar_ft.csv --column ${BUILD_ATMOSPHERE_DATA}/string_ft.csv --missing 9999 --humidity 62.5 --levels 0,1,3,10,20 --z0_from integrated/fit_z0_known.json --top 5000 -o integrated/build_atmosphere_z0_from.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_z0_from PROPERTIES FIXTURES_SETUP build_atmosphere_z0_from_runs FIXTURES_REQUIRED fit_z0_runs)
add_test(NAME build_atmosphere_z0_from_is_harness_profile COMMAND golden_compare integrated/build_atmosphere_z0_from.csv tests/data/atmosphere_builder_expected_ft.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_z0_from_is_harness_profile PROPERTIES FIXTURES_REQUIRED build_atmosphere_z0_from_runs)
add_test(NAME build_atmosphere_z0_twice_refused COMMAND niceops_build_atmosphere ${BUILD_ATMOSPHERE_HARNESS} --z0_from integrated/fit_z0_known.json -o integrated/build_atmosphere_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_z0_twice_refused PROPERTIES PASS_REGULAR_EXPRESSION "give --z0 or --z0_from, not both")
# A fit that reached an end of its search is a limit, and building from it says so; a fit
# that would take no runs is refused.
add_test(NAME build_atmosphere_z0_from_at_bound_warns COMMAND niceops_build_atmosphere build --column ${BUILD_ATMOSPHERE_DATA}/lidar_ft.csv --column ${BUILD_ATMOSPHERE_DATA}/string_ft.csv --missing 9999 --humidity 62.5 --levels 0,1,3,10,20 --z0_from ${BUILD_ATMOSPHERE_DATA}/fit_at_bound.json --top 5000 -o integrated/build_atmosphere_at_bound.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_z0_from_at_bound_warns PROPERTIES PASS_REGULAR_EXPRESSION "that fit came to an end of its search")
add_test(NAME fit_z0_needs_a_run COMMAND niceops_build_atmosphere fit_z0 ${BUILD_ATMOSPHERE_DATA}/fit_run6.csv --missing 9999 --min_runs 0 -o integrated/fit_z0_none.json WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(fit_z0_needs_a_run PROPERTIES PASS_REGULAR_EXPRESSION "--min_runs: Value 0 not in range")
# What the builder refuses, by message.
add_test(NAME build_atmosphere_z0_missing_refused COMMAND niceops_build_atmosphere build --column ${BUILD_ATMOSPHERE_DATA}/lidar_ft.csv --humidity 50 --temperature 288 --top 5000 -o integrated/build_atmosphere_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_z0_missing_refused PROPERTIES PASS_REGULAR_EXPRESSION "give the roughness length, --z0 or --z0_from")
add_test(NAME build_atmosphere_top_missing_refused COMMAND niceops_build_atmosphere build --column ${BUILD_ATMOSPHERE_DATA}/lidar_ft.csv --humidity 50 --temperature 288 --z0 0.1 -o integrated/build_atmosphere_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_top_missing_refused PROPERTIES PASS_REGULAR_EXPRESSION "give --top")
add_test(NAME build_atmosphere_bare_speed_refused COMMAND niceops_build_atmosphere build --wind 10,270 --humidity 50 --temperature 288 --z0 0.1 --top 3000 -o integrated/build_atmosphere_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_bare_speed_refused PROPERTIES PASS_REGULAR_EXPRESSION "--wind: name the speed's unit")
add_test(NAME build_atmosphere_calm_gate_refused COMMAND niceops_build_atmosphere build --column ${BUILD_ATMOSPHERE_DATA}/calm_gate_ft.csv --humidity 50 --temperature 288 --z0 0.1 --top 3000 -o integrated/build_atmosphere_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_calm_gate_refused PROPERTIES PASS_REGULAR_EXPRESSION "gives a calm wind at 125 ft, above the ground")
add_test(NAME build_atmosphere_z0_above_gate_refused COMMAND niceops_build_atmosphere build --wind 10kt,270,5 --humidity 50 --temperature 288 --z0 8 --top 5000 -o integrated/build_atmosphere_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_z0_above_gate_refused PROPERTIES PASS_REGULAR_EXPRESSION "the lowest wind level \\(5 ft\\) has to stand above the roughness length \\(8 ft\\)")
add_test(NAME build_atmosphere_pressure_rises_refused COMMAND niceops_build_atmosphere build --column ${BUILD_ATMOSPHERE_DATA}/pressure_rises_ft.csv --z0 0.1 --top 3000 -o integrated/build_atmosphere_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_pressure_rises_refused PROPERTIES PASS_REGULAR_EXPRESSION "pressure from .*pressure_rises_ft.csv rises with height")
add_test(NAME build_atmosphere_duplicate_source_refused COMMAND niceops_build_atmosphere build --column ${BUILD_ATMOSPHERE_DATA}/string_ft.csv --column ${BUILD_ATMOSPHERE_DATA}/pressure_rises_ft.csv --z0 0.1 --top 3000 -o integrated/build_atmosphere_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_duplicate_source_refused PROPERTIES PASS_REGULAR_EXPRESSION "the temperature comes from both .* name the one to take with --temperature_from")
add_test(NAME build_atmosphere_depth_with_gates_refused COMMAND niceops_build_atmosphere ${BUILD_ATMOSPHERE_HARNESS} --surface_layer_depth 3000 -o integrated/build_atmosphere_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_depth_with_gates_refused PROPERTIES PASS_REGULAR_EXPRESSION "--surface_layer_depth carries the law above a single wind level")
add_test(NAME build_atmosphere_humidity_over_100_refused COMMAND niceops_build_atmosphere build --column ${BUILD_ATMOSPHERE_DATA}/humid_ft.csv --z0 0.1 --top 3000 -o integrated/build_atmosphere_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_humidity_over_100_refused PROPERTIES PASS_REGULAR_EXPRESSION "relative humidity of 103 percent at 100 ft, above saturation")
add_test(NAME build_atmosphere_fill_value_refused COMMAND niceops_build_atmosphere build --column ${BUILD_ATMOSPHERE_DATA}/lidar_ft.csv --humidity 50 --temperature 288 --z0 0.1 --top 5000 -o integrated/build_atmosphere_refused.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(build_atmosphere_fill_value_refused PROPERTIES PASS_REGULAR_EXPRESSION "data row 5 gives a wind speed of 5143.93 m/s, which is not plausible. List a fill value with --missing")
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
set(TERRAIN_BASE -d databases/B407_fixed_load_factor.nod -v databases/B407.json --humidity 20 -t tracks/S-76D/197421AC.csv --origin 39.81043 -74.42631 --frame_bearing 210)
set(TERRAIN_DEMO ${TERRAIN_BASE} -L warn)
set(TERRAIN_PLATE --ground_effect --impedance_model variable_porosity --flow_resistance 200 --plate_table plates/gr1425_vp200.nc)
set(TERRAIN_FLAT --terrain integrated/terrain_flat.nc --terrain_offset -1250 --ground_elevation 54.2651)
set(TERRAIN_HILLS --terrain integrated/terrain_hills.nc --ground_elevation 54.2651)
add_test(NAME terrain_files COMMAND niceops_unit_tests "[.terrain_files]" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_files PROPERTIES FIXTURES_SETUP terrain_file)
add_test(NAME terrain_flat_point_reference COMMAND niceops ${TERRAIN_DEMO} --ground_effect --gridstep 1968 -o integrated/terrain_flat_point_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_point COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_FLAT} --ground_effect --gridstep 1968 -o integrated/terrain_flat_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_plate_reference COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_PLATE} --gridstep 1968 -o integrated/terrain_flat_plate_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_plate COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_FLAT} ${TERRAIN_PLATE} --gridstep 1968 -o integrated/terrain_flat_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_scattered_reference COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_PLATE} -b tests/data/coyle_approach_level.csv -o integrated/terrain_flat_scattered_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_scattered COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_FLAT} ${TERRAIN_PLATE} -b arrays/coyle_approach.csv --observers_on_terrain -o integrated/terrain_flat_scattered.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_flat_point_reference terrain_flat_plate_reference terrain_flat_scattered_reference PROPERTIES FIXTURES_SETUP terrain_flat_runs)
set_tests_properties(terrain_flat_point terrain_flat_plate terrain_flat_scattered PROPERTIES FIXTURES_REQUIRED terrain_file FIXTURES_SETUP terrain_flat_runs)
add_test(NAME terrain_flat_point_is_flat COMMAND golden_compare integrated/terrain_flat_point.csv integrated/terrain_flat_point_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_plate_is_flat COMMAND golden_compare integrated/terrain_flat_plate.csv integrated/terrain_flat_plate_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_scattered_is_flat COMMAND golden_compare integrated/terrain_flat_scattered.csv integrated/terrain_flat_scattered_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_flat_point_is_flat terrain_flat_plate_is_flat terrain_flat_scattered_is_flat PROPERTIES FIXTURES_REQUIRED terrain_flat_runs)
add_test(NAME terrain_hills_point COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} --ground_effect --gridstep 492 -o integrated/terrain_hills_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
# Paths behind the hills are not screened, and the run says so as a warning; over flat terrain it
# has nothing to say.
add_test(NAME terrain_hills_warns_unscreened COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} --ground_effect --gridstep 4000 -o integrated/terrain_hills_warns.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_hills_warns_unscreened PROPERTIES FIXTURES_REQUIRED terrain_file PASS_REGULAR_EXPRESSION "observers, paths had a source or receiver below their mean ground plane")
set_tests_properties(terrain_flat_point PROPERTIES FAIL_REGULAR_EXPRESSION "below their mean ground plane")
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
# Ground classes over the hills (each point receiver's model at its path's own height, over its
# class): one class repeating the command line is the run without classes, byte for byte, and
# a grass and concrete array's peak search equals the march.
set(TERRAIN_CLASSES ${TERRAIN_DEMO} ${TERRAIN_HILLS} --ground_effect --impedance_model variable_porosity --flow_resistance 200 --observers_on_terrain)
add_test(NAME terrain_class_scalar COMMAND niceops ${TERRAIN_CLASSES} -b tests/data/coyle_heights.csv -o integrated/terrain_class_scalar.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_class_single COMMAND niceops ${TERRAIN_CLASSES} --ground_classes tests/data/ground_classes_one.csv -b tests/data/coyle_heights_class7.csv -o integrated/terrain_class_single.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_classes COMMAND niceops ${TERRAIN_CLASSES} --ground_classes tests/data/ground_classes_two.csv -b tests/data/coyle_heights_classes.csv -o integrated/terrain_classes.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_classes_marched COMMAND niceops ${TERRAIN_CLASSES} --ground_classes tests/data/ground_classes_two.csv -b tests/data/coyle_heights_classes.csv -D integrated/terrain_classes_dnl.csv -o integrated/terrain_classes_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_class_scalar terrain_class_single terrain_classes terrain_classes_marched PROPERTIES FIXTURES_REQUIRED terrain_file FIXTURES_SETUP terrain_class_runs)
add_test(NAME ground_class_single_is_scalar_terrain COMMAND ${CMAKE_COMMAND} -E compare_files integrated/terrain_class_single.csv integrated/terrain_class_scalar.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME ground_classes_terrain_windowed_is_exact COMMAND golden_compare integrated/terrain_classes.csv integrated/terrain_classes_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_class_single_is_scalar_terrain ground_classes_terrain_windowed_is_exact PROPERTIES FIXTURES_REQUIRED terrain_class_runs)
# The array's first concrete microphone is on a pole (value 5): it must differ from the all-grass run
# there first, which an admittance taken from the wrong class along a terrain path would not.
add_test(NAME ground_classes_terrain_impedance_is_the_observers COMMAND golden_compare integrated/terrain_classes.csv integrated/terrain_class_scalar.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(ground_classes_terrain_impedance_is_the_observers PROPERTIES FIXTURES_REQUIRED terrain_class_runs PASS_REGULAR_EXPRESSION "mismatch at value 5:")
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
# --propagation_table over --terrain: each pair reads the table at its heights above, and range
# along, its path's mean ground plane, in the frame the ground model takes there.  On level terrain
# that is the flat read, byte for byte (point receivers, plates, and scattered observers stood on
# it); over the hills the peak search must still equal the march, with the ground-model table
# below (rows 0 to 25 ft, Amedee's variable-porosity ground) and with the bump table on plates.
set(TERRAIN_TABLE_PLATE --ground_effect --flow_resistance 225 --plate_table plates/gr1425_db225.nc --propagation_table tables/ground_model_test.nc)
set(TERRAIN_TABLE_GROUND --ground_effect --impedance_model variable_porosity --flow_resistance 200 --propagation_table tables/ground_model_terrain_test.nc)
foreach(height RANGE 0 1400 20)
  list(APPEND TERRAIN_TABLE_SOURCES ${height})
endforeach()
add_test(NAME terrain_table_build COMMAND niceops_propagation_table --method ground_model --database databases/B407_fixed_load_factor.nod --sound_speed 335.636 --impedance_model variable_porosity --flow_resistance 200
    --source_heights ${TERRAIN_TABLE_SOURCES} --receiver_heights 0 3.93701 8 16 25 --max_range 15000 --step 20 -o tables/ground_model_terrain_test.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_table_build PROPERTIES FIXTURES_SETUP terrain_table_file)
add_test(NAME terrain_flat_table_point_reference COMMAND niceops ${TERRAIN_DEMO} --propagation_table tables/ground_model_test.nc --gridstep 1968 -o integrated/terrain_flat_table_point_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_table_point COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_FLAT} --propagation_table tables/ground_model_test.nc --gridstep 1968 -o integrated/terrain_flat_table_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_table_plate_reference COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_TABLE_PLATE} --gridstep 1968 -o integrated/terrain_flat_table_plate_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_table_plate COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_FLAT} ${TERRAIN_TABLE_PLATE} --gridstep 1968 -o integrated/terrain_flat_table_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_table_scattered_reference COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_TABLE_PLATE} -b tests/data/coyle_approach_level.csv -o integrated/terrain_flat_table_scattered_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_table_scattered COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_FLAT} ${TERRAIN_TABLE_PLATE} -b arrays/coyle_approach.csv --observers_on_terrain -o integrated/terrain_flat_table_scattered.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_flat_table_point_reference terrain_flat_table_plate_reference terrain_flat_table_scattered_reference PROPERTIES FIXTURES_REQUIRED propagation_table_file FIXTURES_SETUP terrain_flat_table_runs)
set_tests_properties(terrain_flat_table_point terrain_flat_table_plate terrain_flat_table_scattered PROPERTIES FIXTURES_REQUIRED "terrain_file;propagation_table_file" FIXTURES_SETUP terrain_flat_table_runs FAIL_REGULAR_EXPRESSION "below their mean ground plane")
add_test(NAME terrain_flat_table_point_is_flat COMMAND ${CMAKE_COMMAND} -E compare_files integrated/terrain_flat_table_point.csv integrated/terrain_flat_table_point_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_table_plate_is_flat COMMAND ${CMAKE_COMMAND} -E compare_files integrated/terrain_flat_table_plate.csv integrated/terrain_flat_table_plate_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_table_scattered_is_flat COMMAND ${CMAKE_COMMAND} -E compare_files integrated/terrain_flat_table_scattered.csv integrated/terrain_flat_table_scattered_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_flat_table_point_is_flat terrain_flat_table_plate_is_flat terrain_flat_table_scattered_is_flat PROPERTIES FIXTURES_REQUIRED terrain_flat_table_runs)
add_test(NAME terrain_hills_table_point COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} ${TERRAIN_TABLE_GROUND} --gridstep 492 -o integrated/terrain_hills_table_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_table_point_marched COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} ${TERRAIN_TABLE_GROUND} --gridstep 492 -D integrated/terrain_hills_table_point_dnl.csv -o integrated/terrain_hills_table_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_table_plate COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} ${TERRAIN_TABLE_GROUND} --plate_table plates/gr1425_vp200.nc --gridstep 492 -o integrated/terrain_hills_table_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_table_plate_marched COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} ${TERRAIN_TABLE_GROUND} --plate_table plates/gr1425_vp200.nc --gridstep 492 -D integrated/terrain_hills_table_plate_dnl.csv -o integrated/terrain_hills_table_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_hills_table_point terrain_hills_table_point_marched terrain_hills_table_plate terrain_hills_table_plate_marched PROPERTIES FIXTURES_REQUIRED "terrain_file;terrain_table_file" FIXTURES_SETUP terrain_hills_table_runs)
add_test(NAME terrain_hills_table_bump_plate COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} --ground_effect --flow_resistance 225 --plate_table plates/gr1425_db225.nc --propagation_table tables/bump_test.nc --gridstep 492 -o integrated/terrain_hills_table_bump_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_table_bump_plate_marched COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} --ground_effect --flow_resistance 225 --plate_table plates/gr1425_db225.nc --propagation_table tables/bump_test.nc --gridstep 492 -D integrated/terrain_hills_table_bump_plate_dnl.csv -o integrated/terrain_hills_table_bump_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_hills_table_bump_plate terrain_hills_table_bump_plate_marched PROPERTIES FIXTURES_REQUIRED "terrain_file;propagation_table_bump_file" FIXTURES_SETUP terrain_hills_table_runs)
add_test(NAME terrain_hills_table_point_windowed_is_exact COMMAND golden_compare integrated/terrain_hills_table_point.csv integrated/terrain_hills_table_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_table_plate_windowed_is_exact COMMAND golden_compare integrated/terrain_hills_table_plate.csv integrated/terrain_hills_table_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_table_bump_plate_windowed_is_exact COMMAND golden_compare integrated/terrain_hills_table_bump_plate.csv integrated/terrain_hills_table_bump_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_hills_table_point_windowed_is_exact terrain_hills_table_plate_windowed_is_exact terrain_hills_table_bump_plate_windowed_is_exact PROPERTIES FIXTURES_REQUIRED terrain_hills_table_runs)
# A point receiver's height above a hillside's mean plane is not its pole's: a table with one
# receiver row is warned about up front, and the pairs whose height falls outside the table's rows
# take the ground model, counted and named, never clamped onto a row.
add_test(NAME terrain_table_warns_one_row COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} --propagation_table tables/ground_model_test.nc --gridstep 4000 -o integrated/terrain_table_rows.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_table_warns_one_row PROPERTIES FIXTURES_REQUIRED "terrain_file;propagation_table_file" PASS_REGULAR_EXPRESSION "but it has one receiver row")
add_test(NAME terrain_table_warns_receiver_rows COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} --propagation_table tables/ground_model_test.nc --gridstep 4000 -o integrated/terrain_table_receiver_rows.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_table_warns_receiver_rows PROPERTIES FIXTURES_REQUIRED "terrain_file;propagation_table_file" PASS_REGULAR_EXPRESSION "point-receiver evaluations put the microphone between .* above their path's mean ground plane, outside the table's receiver heights")
# Pairs the table reads over the plane count toward the raised-height warning, with no ground
# model at all: the table alone takes them onto the plane.
add_test(NAME terrain_hills_table_warns_unscreened COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} --no-ground_effect --propagation_table tables/ground_model_terrain_test.nc --receiver_height 3.93701 --gridstep 4000 -o integrated/terrain_hills_table_warns.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_hills_table_warns_unscreened PROPERTIES FIXTURES_REQUIRED "terrain_file;terrain_table_file" PASS_REGULAR_EXPRESSION "observers, paths had a source or receiver below their mean ground plane")
# ... and only those: a pair the table misses, with no ground model to fall back on, is free field
# and never used the plane.  Every point receiver here misses the one-row table, and counting
# those pairs as well raised the same warning at 3 of 9 observers.
add_test(NAME terrain_hills_table_misses_are_free_field COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} --no-ground_effect --propagation_table tables/ground_model_test.nc --gridstep 4000 -o integrated/terrain_hills_table_misses.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_hills_table_misses_are_free_field PROPERTIES FIXTURES_REQUIRED "terrain_file;propagation_table_file" PASS_REGULAR_EXPRESSION "point-receiver evaluations put the microphone .* and took free-field propagation" FAIL_REGULAR_EXPRESSION "below their mean ground plane")
add_test(NAME terrain_table_refuses_refraction COMMAND niceops -d databases/B407_fixed_load_factor.nod -v databases/B407.json -t tracks/S-76D/197421AC.csv --origin 39.81043 -74.42631 --frame_bearing 210 ${TERRAIN_FLAT} --atmosphere atmospheres/B407_uniform.csv --refraction linear --propagation_table tables/ground_model_test.nc -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_table_refuses_refraction PROPERTIES FIXTURES_REQUIRED "terrain_file;propagation_table_file" PASS_REGULAR_EXPRESSION "does not combine with --refraction")
# A table that refracts is read over terrain as if its stratification followed each path's mean
# plane; past 5 degrees of tilt (the 5.7 degree plane) the run warns.  A ground-model table has no
# stratification to tilt and never does.  Plates, which always read the table, carry it.
add_test(NAME terrain_refracting_table_build COMMAND niceops_unit_tests "[.refracting_table]" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_refracting_table_build PROPERTIES FIXTURES_REQUIRED propagation_table_file FIXTURES_SETUP terrain_refracting_table_file)
add_test(NAME terrain_steep_refracting_table_warns COMMAND niceops ${TERRAIN_DEMO} --terrain integrated/terrain_steep.nc --ground_elevation 54.2651 --ground_effect --flow_resistance 225 --plate_table plates/gr1425_db225.nc --propagation_table tables/refracting_test.nc --gridstep 4000 -o integrated/terrain_steep_refracting.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_steep_refracting_table_warns PROPERTIES FIXTURES_REQUIRED "terrain_file;terrain_refracting_table_file" PASS_REGULAR_EXPRESSION "refracts, and over --terrain it is read by height above each path's mean ground plane.*had a plane tilted more than 5 degrees")
add_test(NAME terrain_steep_ground_model_table_is_quiet COMMAND niceops ${TERRAIN_DEMO} --terrain integrated/terrain_steep.nc --ground_elevation 54.2651 ${TERRAIN_TABLE_PLATE} --gridstep 4000 -o integrated/terrain_steep_ground_model.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_steep_ground_model_table_is_quiet PROPERTIES FIXTURES_REQUIRED "terrain_file;propagation_table_file" FAIL_REGULAR_EXPRESSION "refracts, and over --terrain")
# ... and a ground-model table over the hills is the terrain ground model, to the table's
# interpolation, which is the plumbing check: it fails if the read takes the source's height
# above z = 0 for its height above the plane.  A plate reads the five-row table's row at 0 plus K,
# the plate in the plane's frame (0.014 dB at worst; FIDELITY, Propagation tables).  A point
# receiver reads between rows at its own height above the plane, where the ground's fringes in
# receiver height are finer than any practical row spacing above a few hundred hertz, so its check
# runs from 20 to 250 Hz over rows every half foot (0.17 dB at worst).
set(TERRAIN_LOW_TABLE_ROWS 0)
foreach(tenths RANGE 5 220 5)
  math(EXPR feet "${tenths} / 10")
  math(EXPR fraction "${tenths} % 10")
  list(APPEND TERRAIN_LOW_TABLE_ROWS ${feet}.${fraction})
endforeach()
add_test(NAME terrain_low_table_build COMMAND niceops_propagation_table --method ground_model --database databases/B407_fixed_load_factor.nod --band_limits 20 250 --sound_speed 335.636 --impedance_model variable_porosity --flow_resistance 200
    --source_heights ${TERRAIN_TABLE_SOURCES} --receiver_heights ${TERRAIN_LOW_TABLE_ROWS} --max_range 15000 --step 20 -o tables/ground_model_terrain_low_test.nc
    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_low_table_build PROPERTIES FIXTURES_SETUP terrain_low_table_file)
add_test(NAME terrain_hills_low_table_point COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} --ground_effect --impedance_model variable_porosity --flow_resistance 200 --band_limits 20 250 --propagation_table tables/ground_model_terrain_low_test.nc --gridstep 492 -o integrated/terrain_hills_low_table_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_hills_low_table_point PROPERTIES FIXTURES_REQUIRED "terrain_file;terrain_low_table_file" FIXTURES_SETUP terrain_hills_low_runs FAIL_REGULAR_EXPRESSION "outside the table's receiver heights")
add_test(NAME terrain_hills_low_point_reference COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_HILLS} --ground_effect --impedance_model variable_porosity --flow_resistance 200 --band_limits 20 250 --exact_ground --gridstep 492 -o integrated/terrain_hills_low_point_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_hills_low_point_reference PROPERTIES FIXTURES_REQUIRED terrain_file FIXTURES_SETUP terrain_hills_low_runs)
add_test(NAME terrain_hills_table_point_is_the_ground_model COMMAND golden_compare integrated/terrain_hills_low_table_point.csv integrated/terrain_hills_low_point_reference.csv 0.2 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_hills_table_point_is_the_ground_model PROPERTIES FIXTURES_REQUIRED terrain_hills_low_runs)
add_test(NAME terrain_hills_table_plate_is_the_plate_model COMMAND golden_compare integrated/terrain_hills_table_plate.csv integrated/terrain_hills_plate.csv 0.02 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_hills_table_plate_is_the_plate_model PROPERTIES FIXTURES_REQUIRED "terrain_hills_table_runs;terrain_hills_runs")
# A plate whose source's foot lies behind it on its path's mean line (d' < 0, nearly overhead on
# ground falling toward the source) takes K at the turned separation's azimuth.  Those pairs are
# too few in the grids above to move them: K from the true separation's azimuth passed every other
# test, and moves these plates, under the approach on the 5.7 degree plane, by 2e-4 to 7e-4 dB.
add_test(NAME terrain_steep_table_plate COMMAND niceops ${TERRAIN_DEMO} --terrain integrated/terrain_steep.nc --ground_elevation 54.2651 ${TERRAIN_TABLE_PLATE} -b tests/data/steep_under_track.csv --observers_on_terrain -o integrated/terrain_steep_table_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_steep_table_plate PROPERTIES FIXTURES_REQUIRED "terrain_file;propagation_table_file" FIXTURES_SETUP terrain_steep_table_plate_run)
add_test(NAME terrain_steep_table_plate_golden COMMAND golden_compare integrated/terrain_steep_table_plate.csv golden/terrain_steep_table_plate.csv 1e-5 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_steep_table_plate_golden PROPERTIES FIXTURES_REQUIRED terrain_steep_table_plate_run)
# --terrain_screening cnossos (src/screening.hpp; THEORY 6.4).  A run whose info log the test
# matches on goes through tests/run_logged.cmake, which gives it its own log and prints it.  A
# PASS_REGULAR_EXPRESSION makes CTest ignore the exit status, so a run that exits nonzero fails on
# the line run_logged prints for it: a FAIL_REGULAR_EXPRESSION set later must keep it.
set(LOGGED_RUN_FAILED "run_logged: the command exited with")
function(add_logged_run name)
  add_test(NAME ${name} COMMAND ${CMAKE_COMMAND} -DLOG=integrated/${name}.log
      "-DCOMMAND=$<TARGET_FILE:niceops>;${ARGN};-L;info;-l;integrated/${name}.log" -P tests/run_logged.cmake
      WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(${name} PROPERTIES FAIL_REGULAR_EXPRESSION "${LOGGED_RUN_FAILED}")
endfunction()
set(SCREENING --terrain_screening cnossos)
set(TERRAIN_RIDGE --terrain integrated/terrain_ridge.nc --terrain_offset -1250 --ground_elevation 54.2651)
set(TERRAIN_TILTED --terrain integrated/terrain_tilted.nc --terrain_offset -1250 --ground_elevation 54.2651)
# Level terrain screened is the run without terrain, byte for byte, with no pair screened; so is a
# plane falling 4% under the approach, against itself unscreened; and so are the hills, where no
# sample rises above any path's line of sight (the hills run is terrain_hills_point's).
add_logged_run(terrain_flat_point_screened ${TERRAIN_BASE} ${TERRAIN_FLAT} --ground_effect ${SCREENING} --gridstep 1968 -o integrated/terrain_flat_point_screened.csv)
add_logged_run(terrain_flat_plate_screened ${TERRAIN_BASE} ${TERRAIN_FLAT} ${TERRAIN_PLATE} ${SCREENING} --gridstep 1968 -o integrated/terrain_flat_plate_screened.csv)
add_logged_run(terrain_flat_scattered_screened ${TERRAIN_BASE} ${TERRAIN_FLAT} ${TERRAIN_PLATE} ${SCREENING} -b arrays/coyle_approach.csv --observers_on_terrain -o integrated/terrain_flat_scattered_screened.csv)
add_test(NAME terrain_tilted_point COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_TILTED} --ground_effect --gridstep 1968 -o integrated/terrain_tilted_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_logged_run(terrain_tilted_point_screened ${TERRAIN_BASE} ${TERRAIN_TILTED} --ground_effect ${SCREENING} --gridstep 1968 -o integrated/terrain_tilted_point_screened.csv)
add_test(NAME terrain_tilted_plate COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_TILTED} ${TERRAIN_PLATE} --gridstep 1968 -o integrated/terrain_tilted_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_logged_run(terrain_tilted_plate_screened ${TERRAIN_BASE} ${TERRAIN_TILTED} ${TERRAIN_PLATE} ${SCREENING} --gridstep 1968 -o integrated/terrain_tilted_plate_screened.csv)
add_logged_run(terrain_hills_point_screened ${TERRAIN_BASE} ${TERRAIN_HILLS} --ground_effect ${SCREENING} --gridstep 492 -o integrated/terrain_hills_point_screened.csv)
set_tests_properties(terrain_flat_point_screened terrain_flat_plate_screened terrain_flat_scattered_screened
    terrain_tilted_point terrain_tilted_point_screened terrain_tilted_plate terrain_tilted_plate_screened terrain_hills_point_screened
    PROPERTIES FIXTURES_REQUIRED terrain_file FIXTURES_SETUP terrain_unscreened_runs)
set_tests_properties(terrain_flat_point_screened terrain_flat_plate_screened terrain_flat_scattered_screened
    terrain_tilted_point_screened terrain_tilted_plate_screened terrain_hills_point_screened
    PROPERTIES PASS_REGULAR_EXPRESSION "--terrain_screening: 0 evaluations were screened"
    FAIL_REGULAR_EXPRESSION "--terrain_screening: [1-9];${LOGGED_RUN_FAILED}")
add_test(NAME terrain_flat_point_screened_is_flat COMMAND ${CMAKE_COMMAND} -E compare_files integrated/terrain_flat_point_screened.csv integrated/terrain_flat_point_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_plate_screened_is_flat COMMAND ${CMAKE_COMMAND} -E compare_files integrated/terrain_flat_plate_screened.csv integrated/terrain_flat_plate_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_flat_scattered_screened_is_flat COMMAND ${CMAKE_COMMAND} -E compare_files integrated/terrain_flat_scattered_screened.csv integrated/terrain_flat_scattered_reference.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_tilted_point_is_unscreened COMMAND ${CMAKE_COMMAND} -E compare_files integrated/terrain_tilted_point_screened.csv integrated/terrain_tilted_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_tilted_plate_is_unscreened COMMAND ${CMAKE_COMMAND} -E compare_files integrated/terrain_tilted_plate_screened.csv integrated/terrain_tilted_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_hills_unblocked_is_unscreened COMMAND ${CMAKE_COMMAND} -E compare_files integrated/terrain_hills_point_screened.csv integrated/terrain_hills_point.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_flat_point_screened_is_flat terrain_flat_plate_screened_is_flat terrain_flat_scattered_screened_is_flat
    terrain_tilted_point_is_unscreened terrain_tilted_plate_is_unscreened terrain_hills_unblocked_is_unscreened
    PROPERTIES FIXTURES_REQUIRED "terrain_unscreened_runs;terrain_flat_runs;terrain_hills_runs")
# Behind the ridge the peak search, with the screening's bound term, still equals the march: point
# receivers, plates, no ground model, mixed ground classes, and a propagation table, which a
# blocked pair leaves for the ground model (logged).
add_logged_run(terrain_ridge_screened_point ${TERRAIN_BASE} ${TERRAIN_RIDGE} --ground_effect ${SCREENING} --gridstep 492 -o integrated/terrain_ridge_screened_point.csv)
add_test(NAME terrain_ridge_screened_point_marched COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_RIDGE} --ground_effect ${SCREENING} --gridstep 492 -D integrated/terrain_ridge_screened_point_dnl.csv -o integrated/terrain_ridge_screened_point_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_logged_run(terrain_ridge_screened_plate ${TERRAIN_BASE} ${TERRAIN_RIDGE} ${TERRAIN_PLATE} ${SCREENING} --gridstep 492 -o integrated/terrain_ridge_screened_plate.csv)
add_test(NAME terrain_ridge_screened_plate_marched COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_RIDGE} ${TERRAIN_PLATE} ${SCREENING} --gridstep 492 -D integrated/terrain_ridge_screened_plate_dnl.csv -o integrated/terrain_ridge_screened_plate_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_logged_run(terrain_ridge_screened_noground ${TERRAIN_BASE} ${TERRAIN_RIDGE} --no-ground_effect ${SCREENING} --gridstep 492 -o integrated/terrain_ridge_screened_noground.csv)
add_test(NAME terrain_ridge_screened_noground_marched COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_RIDGE} --no-ground_effect ${SCREENING} --gridstep 492 -D integrated/terrain_ridge_screened_noground_dnl.csv -o integrated/terrain_ridge_screened_noground_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_logged_run(terrain_ridge_screened_classes ${TERRAIN_BASE} ${TERRAIN_RIDGE} --ground_effect --impedance_model variable_porosity --flow_resistance 200 --observers_on_terrain ${SCREENING} --ground_classes tests/data/ground_classes_two.csv -b tests/data/coyle_heights_classes.csv -o integrated/terrain_ridge_screened_classes.csv)
add_test(NAME terrain_ridge_screened_classes_marched COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_RIDGE} --ground_effect --impedance_model variable_porosity --flow_resistance 200 --observers_on_terrain ${SCREENING} --ground_classes tests/data/ground_classes_two.csv -b tests/data/coyle_heights_classes.csv -D integrated/terrain_ridge_screened_classes_dnl.csv -o integrated/terrain_ridge_screened_classes_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_ridge_screened_point terrain_ridge_screened_plate terrain_ridge_screened_noground terrain_ridge_screened_classes
    PROPERTIES PASS_REGULAR_EXPRESSION "--terrain_screening: [1-9][0-9]* evaluations were screened")
set_tests_properties(terrain_ridge_screened_point terrain_ridge_screened_point_marched terrain_ridge_screened_plate terrain_ridge_screened_plate_marched
    terrain_ridge_screened_noground terrain_ridge_screened_noground_marched terrain_ridge_screened_classes terrain_ridge_screened_classes_marched
    PROPERTIES FIXTURES_REQUIRED terrain_file FIXTURES_SETUP terrain_ridge_runs)
add_logged_run(terrain_ridge_screened_table ${TERRAIN_BASE} ${TERRAIN_RIDGE} ${TERRAIN_TABLE_GROUND} ${SCREENING} --gridstep 492 -o integrated/terrain_ridge_screened_table.csv)
add_test(NAME terrain_ridge_screened_table_marched COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_RIDGE} ${TERRAIN_TABLE_GROUND} ${SCREENING} --gridstep 492 -D integrated/terrain_ridge_screened_table_dnl.csv -o integrated/terrain_ridge_screened_table_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_ridge_screened_table PROPERTIES PASS_REGULAR_EXPRESSION " [1-9][0-9]* of the screened evaluations were not read from the propagation table")
set_tests_properties(terrain_ridge_screened_table terrain_ridge_screened_table_marched PROPERTIES FIXTURES_REQUIRED "terrain_file;terrain_table_file" FIXTURES_SETUP terrain_ridge_runs)
foreach(kind point plate noground classes table)
  add_test(NAME terrain_ridge_screened_${kind}_windowed_is_exact COMMAND golden_compare integrated/terrain_ridge_screened_${kind}.csv integrated/terrain_ridge_screened_${kind}_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
  set_tests_properties(terrain_ridge_screened_${kind}_windowed_is_exact PROPERTIES FIXTURES_REQUIRED terrain_ridge_runs)
endforeach()
# The windowed/marched pairs pass with the screening term left out of the bound, so a coarse
# ridge run checks, from its log, that the term raises the bound in some band.
# An observer file's own z, kept without --observers_on_terrain, is where the
# microphone is evaluated, and so where its horizon is taken from: the same
# observers 3000 ft above their surveyed ground see over the ridge, which on
# the terrain screens 9 of them.
add_logged_run(terrain_ridge_screened_file_z ${TERRAIN_BASE} ${TERRAIN_RIDGE} --ground_effect ${SCREENING} -b tests/data/coyle_heights_aloft.csv -o integrated/terrain_ridge_screened_file_z.csv)
set_tests_properties(terrain_ridge_screened_file_z PROPERTIES FIXTURES_REQUIRED terrain_file
    PASS_REGULAR_EXPRESSION "--terrain_screening: 0 evaluations were screened")
add_logged_run(terrain_ridge_screened_bound ${TERRAIN_BASE} ${TERRAIN_RIDGE} --ground_effect ${SCREENING} --gridstep 1968 -o integrated/terrain_ridge_screened_bound.csv)
set_tests_properties(terrain_ridge_screened_bound PROPERTIES FIXTURES_REQUIRED terrain_file
    PASS_REGULAR_EXPRESSION "bound's screened term is above its ground term in [1-9][0-9]* of [0-9]+ bands")
# The windowed/marched pairs compare one binary with itself, so the ridge runs are pinned (made by
# the build that added --terrain_screening; regenerate only with a deliberate change).
add_test(NAME terrain_ridge_screened_point_golden COMMAND golden_compare integrated/terrain_ridge_screened_point.csv golden/terrain_ridge_screened_point.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_ridge_screened_plate_golden COMMAND golden_compare integrated/terrain_ridge_screened_plate.csv golden/terrain_ridge_screened_plate.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_ridge_screened_point_golden terrain_ridge_screened_plate_golden PROPERTIES FIXTURES_REQUIRED terrain_ridge_runs)
# --terrain_screening kirchhoff (src/screening.hpp; THEORY 6.4): a blocked pair keeps its own ground,
# plate or table and loses a knife edge's Fresnel-Kirchhoff loss, at least 10 lg 4.  Where nothing is
# blocked it is the unscreened run byte for byte, shown on the falling plane.  It shares CNOSSOS's
# gate (planes->blocked, then the rubber band), so level terrain and the hills, which the CNOSSOS
# runs above take through that gate, are not run again here.
set(KIRCHHOFF --terrain_screening kirchhoff)
add_logged_run(terrain_tilted_plate_kirchhoff ${TERRAIN_BASE} ${TERRAIN_TILTED} ${TERRAIN_PLATE} ${KIRCHHOFF} --gridstep 1968 -o integrated/terrain_tilted_plate_kirchhoff.csv)
set_tests_properties(terrain_tilted_plate_kirchhoff
    PROPERTIES FIXTURES_REQUIRED terrain_file FIXTURES_SETUP terrain_kirchhoff_unscreened_runs
    PASS_REGULAR_EXPRESSION "--terrain_screening: 0 evaluations were screened"
    FAIL_REGULAR_EXPRESSION "--terrain_screening: [1-9];${LOGGED_RUN_FAILED}")
add_test(NAME terrain_tilted_plate_kirchhoff_is_unscreened COMMAND ${CMAKE_COMMAND} -E compare_files integrated/terrain_tilted_plate_kirchhoff.csv integrated/terrain_tilted_plate.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_tilted_plate_kirchhoff_is_unscreened
    PROPERTIES FIXTURES_REQUIRED "terrain_kirchhoff_unscreened_runs;terrain_unscreened_runs")
# Behind the ridge the peak search, with no bound term for the knife edge, equals the march.  The
# knife edge adds no bound term (a blocked pair is its unscreened level less a loss), so one pair, the
# cheapest, shows it; the point, plate and table call sites are pinned by their goldens below.
add_logged_run(terrain_ridge_kirchhoff_point ${TERRAIN_BASE} ${TERRAIN_RIDGE} --ground_effect ${KIRCHHOFF} --gridstep 492 -o integrated/terrain_ridge_kirchhoff_point.csv)
add_logged_run(terrain_ridge_kirchhoff_plate ${TERRAIN_BASE} ${TERRAIN_RIDGE} ${TERRAIN_PLATE} ${KIRCHHOFF} --gridstep 492 -o integrated/terrain_ridge_kirchhoff_plate.csv)
add_logged_run(terrain_ridge_kirchhoff_noground ${TERRAIN_BASE} ${TERRAIN_RIDGE} --no-ground_effect ${KIRCHHOFF} --gridstep 492 -o integrated/terrain_ridge_kirchhoff_noground.csv)
add_test(NAME terrain_ridge_kirchhoff_noground_marched COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_RIDGE} --no-ground_effect ${KIRCHHOFF} --gridstep 492 -D integrated/terrain_ridge_kirchhoff_noground_dnl.csv -o integrated/terrain_ridge_kirchhoff_noground_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_logged_run(terrain_ridge_kirchhoff_table ${TERRAIN_BASE} ${TERRAIN_RIDGE} ${TERRAIN_TABLE_GROUND} ${KIRCHHOFF} --gridstep 492 -o integrated/terrain_ridge_kirchhoff_table.csv)
set_tests_properties(terrain_ridge_kirchhoff_point terrain_ridge_kirchhoff_plate terrain_ridge_kirchhoff_noground
    PROPERTIES PASS_REGULAR_EXPRESSION "--terrain_screening: [1-9][0-9]* evaluations were screened")
set_tests_properties(terrain_ridge_kirchhoff_table PROPERTIES PASS_REGULAR_EXPRESSION "The screened evaluations kept the propagation table")
set_tests_properties(terrain_ridge_kirchhoff_point terrain_ridge_kirchhoff_plate
    terrain_ridge_kirchhoff_noground terrain_ridge_kirchhoff_noground_marched
    PROPERTIES FIXTURES_REQUIRED terrain_file FIXTURES_SETUP terrain_ridge_kirchhoff_runs)
set_tests_properties(terrain_ridge_kirchhoff_table PROPERTIES FIXTURES_REQUIRED "terrain_file;terrain_table_file" FIXTURES_SETUP terrain_ridge_kirchhoff_runs)
add_test(NAME terrain_ridge_kirchhoff_noground_windowed_is_exact COMMAND golden_compare integrated/terrain_ridge_kirchhoff_noground.csv integrated/terrain_ridge_kirchhoff_noground_marched.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_ridge_kirchhoff_noground_windowed_is_exact PROPERTIES FIXTURES_REQUIRED terrain_ridge_kirchhoff_runs)
# ... the knife edge acts (the ridge's plates differ from CNOSSOS's: golden_compare's mismatch is the
# pass), and the runs are pinned (made by the build that added kirchhoff; regenerate only with a
# deliberate change).
add_test(NAME terrain_ridge_kirchhoff_is_not_cnossos COMMAND golden_compare integrated/terrain_ridge_kirchhoff_plate.csv integrated/terrain_ridge_screened_plate.csv 0.01 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_ridge_kirchhoff_is_not_cnossos PROPERTIES PASS_REGULAR_EXPRESSION "mismatch at value" FIXTURES_REQUIRED "terrain_ridge_kirchhoff_runs;terrain_ridge_runs")
add_test(NAME terrain_ridge_kirchhoff_point_golden COMMAND golden_compare integrated/terrain_ridge_kirchhoff_point.csv golden/terrain_ridge_kirchhoff_point.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_ridge_kirchhoff_plate_golden COMMAND golden_compare integrated/terrain_ridge_kirchhoff_plate.csv golden/terrain_ridge_kirchhoff_plate.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_ridge_kirchhoff_point_golden terrain_ridge_kirchhoff_plate_golden PROPERTIES FIXTURES_REQUIRED terrain_ridge_kirchhoff_runs)
# The table's and the free field's call sites, each of which a windowed/marched pair passes without
# the loss: pinned too (the table's looser, as it is built by the test run).
add_test(NAME terrain_ridge_kirchhoff_noground_golden COMMAND golden_compare integrated/terrain_ridge_kirchhoff_noground.csv golden/terrain_ridge_kirchhoff_noground.csv 1e-6 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME terrain_ridge_kirchhoff_table_golden COMMAND golden_compare integrated/terrain_ridge_kirchhoff_table.csv golden/terrain_ridge_kirchhoff_table.csv 1e-5 WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_ridge_kirchhoff_noground_golden terrain_ridge_kirchhoff_table_golden PROPERTIES FIXTURES_REQUIRED terrain_ridge_kirchhoff_runs)
add_test(NAME terrain_screening_rejects_unknown COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_RIDGE} --terrain_screening deygout -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_screening_rejects_unknown PROPERTIES PASS_REGULAR_EXPRESSION "--terrain_screening: deygout not in")
# Refusals.
add_test(NAME terrain_screening_needs_terrain COMMAND niceops ${TERRAIN_DEMO} ${SCREENING} -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_screening_needs_terrain PROPERTIES PASS_REGULAR_EXPRESSION "--terrain_screening needs --terrain")
add_test(NAME screening_min_height_needs_screening COMMAND niceops --screening_min_height 1 -o integrated/never.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(screening_min_height_needs_screening PROPERTIES PASS_REGULAR_EXPRESSION "--screening_min_height is for --terrain_screening")
add_test(NAME terrain_screening_needs_spectral COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_RIDGE} --no-spectral --no-ground_effect ${SCREENING} -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_screening_needs_spectral PROPERTIES FIXTURES_REQUIRED terrain_file PASS_REGULAR_EXPRESSION "--terrain_screening is per band, and the broadband path")
add_test(NAME terrain_screening_needs_spectra COMMAND niceops -d databases/S-76D_M3.nod -v databases/S-76D.json -t tracks/S-76D/197421AC.csv --origin 39.81043 -74.42631 --frame_bearing 210 ${TERRAIN_RIDGE} ${SCREENING} --gridstep 4000 -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_screening_needs_spectra PROPERTIES FIXTURES_REQUIRED terrain_file PASS_REGULAR_EXPRESSION "--terrain_screening is per band, and databases/S-76D_M3.nod carries no per-band spectra")
add_test(NAME terrain_screening_table_needs_ground COMMAND niceops ${TERRAIN_DEMO} ${TERRAIN_RIDGE} --no-ground_effect --propagation_table tables/ground_model_test.nc ${SCREENING} -o integrated/never_written.csv WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
set_tests_properties(terrain_screening_table_needs_ground PROPERTIES FIXTURES_REQUIRED "terrain_file;propagation_table_file" PASS_REGULAR_EXPRESSION "--terrain_screening cnossos with --propagation_table needs a ground model")

# Which tests CI runs (label `extended` for the rest), checked against the fixtures.
include(tests/extended.cmake)
