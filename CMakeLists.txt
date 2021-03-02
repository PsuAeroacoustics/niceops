cmake_minimum_required(VERSION 3.1)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-d)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen3/cmake)

include(BuildDExecutable)


project(NICEOPS VERSION 1.0 LANGUAGES C CXX D)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX
    "build" CACHE PATH "default install prefix" FORCE
    )
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

include(ExtractDMDSystemLinker)

#Set dylib paths relative to executable path under OSX
set(CMAKE_MACOSX_RPATH 0)
set(CMAKE_INSTALL_NAME_DIR @executable_path)
#Set attributes for Xcode project
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++14")
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
set(CMAKE_CXX_STANDARD 14)

#include_directories(SYSTEM /usr/local/include /opt/local/include /sw/include /usr/include external/splinter/include external/matplotlibd/source /usr/local/include/eigen3 /usr/include/eigen3 external/numd/source)
include_directories(SYSTEM /usr/local/include /opt/local/include /sw/include /usr/include external/splinter/include external/matplotlibd/source external/eigen3 external/numd/source)
#include_directories(SYSTEM /usr/local/include /opt/local/include /sw/include /usr/include external/splinter/include external/matplotlibd/source external/numd/source)

# Including the source directory is apparently required to support out of source builds of D modules
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_definitions(-D_THREAD_SAFE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
          "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
          FORCE)
endif (NOT CMAKE_BUILD_TYPE)

find_package(NetCDF)
message(NOTICE ${NetCDF_VERSION})

# Force shared version to be used to ensure good linking . (Some systems like Homebrew install static and shared in the same place)
find_file(NCDF NAMES libnetcdf.dylib libnetcdf.so netcdf.lib PATHS /usr/local/lib/ /opt/local/lib/ /sw/lib/ /usr/lib64/ /usr/lib/ /usr/lib/x86_64-linux-gnu/netcdf/mpi "C:\\Program Files\\netCDF ${NetCDF_VERSION}\\lib" "C:\\Program Files (x86)\\netCDF ${NetCDF_VERSION}\\lib")

find_package(Boost)
find_package(CGAL COMPONENTS Core)
include(${CGAL_USE_FILE})

if(NOT ${GMP_INCLUDES} STREQUAL "")
  include_directories(${GMP_INCLUDES})
endif()

if(NOT ${MPFR_INCLUDES} STREQUAL "")
  include_directories(${MPFR_INCLUDES})
endif()

add_library(sphere_database sphere_database.d netcdf.d)
target_link_libraries(sphere_database ${NCDF})
set_target_properties(
  sphere_database PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ./
    LIBRARY_OUTPUT_DIRECTORY ./
    ARCHIVE_OUTPUT_DIRECTORY ./
    ARCHIVE_OUTPUT_NAME sphere_database
    LIBRARY_OUTPUT_NAME sphere_database
    RUNTIME_OUTPUT_NAME sphere_database
    COMPILE_FLAGS ""
    LINKER_LANGUAGE CXX
    LINK_FLAGS "${SANITIZE_LDFLAGS}"
    # Use a custom .props file to set up Visual D (import paths, predefined versions...).
    VS_USER_PROPS "${PROJECT_SOURCE_DIR}/cmake-d/VisualD.props"
)

add_library(interpolator interpolator.cpp)
target_link_libraries(interpolator ${CGAL})

add_subdirectory(external/splinter)
add_library(trajectory trajectory.d wgs84tocartesian.d worldtransforms.d splines.d splinter.cpp)
target_link_libraries(trajectory splinter-static)
set_target_properties(
    trajectory PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ./
    LIBRARY_OUTPUT_DIRECTORY ./
    ARCHIVE_OUTPUT_DIRECTORY ./
    ARCHIVE_OUTPUT_NAME trajectory
    LIBRARY_OUTPUT_NAME trajectory
    RUNTIME_OUTPUT_NAME trajectory
    COMPILE_FLAGS ""
    LINKER_LANGUAGE CXX
    LINK_FLAGS "${SANITIZE_LDFLAGS}"
    # Use a custom .props file to set up Visual D (import paths, predefined versions...).
    VS_USER_PROPS "${PROJECT_SOURCE_DIR}/cmake-d/VisualD.props"
)

add_subdirectory(external/numd)

add_subdirectory(external/matplotlibd)
#function(build_d_executable target_name output_exe d_src_files compiler_args linker_args extra_compile_deps link_deps compile_separately)
#build_d_executable(
#    "${LDMD_EXE}"
#    "${LDMD_EXE_FULL}"
#    "${LDMD_D_SOURCE_FILES}"
#    "${DFLAGS_BUILD_TYPE}"
#    "${LDC_LINKERFLAG_LIST}"
#    ""
#    "LDMD_CXX_LIB"
#    ${COMPILE_D_MODULES_SEPARATELY}
#)
set(views_dir ${CMAKE_CURRENT_SOURCE_DIR}/external/matplotlibd/views)
set(LDC_LINK_MANUALLY TRUE)

build_d_executable(
  niceops
  ./niceops
  ${CMAKE_CURRENT_SOURCE_DIR}/niceops.d
  "-J${views_dir} -Iexternal/matplotlibd/source -Iexternal/numd/source"
  "-Wl,--copy-dt-needed-entries"
  ""
  ""
  "sphere_database;interpolator;trajectory;matplotlib-d;numd"
  TRUE
)
#add_executable(niceops niceops.d)
#target_link_libraries(niceops sphere_database interpolator trajectory matplotlib-d numd)

#set(views_dir ${CMAKE_CURRENT_SOURCE_DIR}/external/matplotlibd/views)
#target_compile_options(niceops PRIVATE "-J ${views_dir}")
