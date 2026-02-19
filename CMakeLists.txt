cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-d)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
#list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen3/cmake)

#include(BuildDExecutable)


project(NICEOPS VERSION 1.0 LANGUAGES C CXX D)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX
    "build" CACHE PATH "default install prefix" FORCE
    )
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

#list(APPEND DFLAGS_BASE "-m64")
#include(ExtractDMDSystemLinker)

if("${CMAKE_CL_64}" STREQUAL "1")
  message(NOTICE "64 bit!")
  set(CMAKE_D_FLAGS "${CMAKE_D_FLAGS}-m64")
  set(CMAKE_D_LINK_FLAGS "${CMAKE_D_LINK_FLAGS}-m64")
  set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS}-m64")
  set(CMAKE_STATIC_LIBRARY_CREATE_D_FLAGS "${CMAKE_STATIC_LIBRARY_CREATE_D_FLAGS}-m64")
  set(CMAKE_SHARED_LIBRARY_D_FLAGS "${CMAKE_SHARED_LIBRARY_D_FLAGS}-m64")
  set(CMAKE_STATIC_LIBRARY_CREATE_D_FLAGS "${CMAKE_STATIC_LIBRARY_CREATE_D_FLAGS}-m64")
  
endif()


#Set dylib paths relative to executable path under OSX
set(CMAKE_MACOSX_RPATH 0)
set(CMAKE_INSTALL_NAME_DIR @executable_path)
#Set attributes for Xcode project
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++17")
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
set(CMAKE_CXX_STANDARD 14)



#include_directories(SYSTEM /usr/local/include /opt/local/include /sw/include /usr/include external/splinter/include external/matplotlibd/source /usr/local/include/eigen3 /usr/include/eigen3 external/numd/source)
include_directories(SYSTEM /usr/local/include /opt/local/include /sw/include /usr/include external/splinter/include external/splinter/thirdparty/Eigen external/matplotlibd/source /usr/local/include/eigen3 /usr/include/eigen3 external/numd/source)
#include_directories(SYSTEM /usr/local/include /opt/local/include /sw/include /usr/include external/splinter/include external/matplotlibd/source external/numd/source)

# Including the source directory is apparently required to support out of source builds of D modules
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_definitions(-D_THREAD_SAFE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
          "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
          FORCE)
endif (NOT CMAKE_BUILD_TYPE)

find_package(NetCDF REQUIRED)
message(NOTICE ${NetCDF_VERSION})

find_package(Boost)
find_package(CGAL COMPONENTS Core REQUIRED)

if(CGAL_USE_FILE)
  include("${CGAL_USE_FILE}")
endif()

if(NOT "${GMP_INCLUDES}" STREQUAL "")
  include_directories(${GMP_INCLUDES})
endif()

if(NOT "${MPFR_INCLUDES}" STREQUAL "")
  include_directories(${MPFR_INCLUDES})
endif()

add_library(sphere_database sphere_database.d netcdf.d)
target_link_libraries(sphere_database NetCDF::NetCDF)
set_target_properties(
  sphere_database PROPERTIES
    #RUNTIME_OUTPUT_DIRECTORY ./
    #LIBRARY_OUTPUT_DIRECTORY ./
    #ARCHIVE_OUTPUT_DIRECTORY ./
    #ARCHIVE_OUTPUT_NAME sphere_database
    #LIBRARY_OUTPUT_NAME sphere_database
    #RUNTIME_OUTPUT_NAME sphere_database
    #COMPILE_FLAGS ""
    LINKER_LANGUAGE CXX
    #LINK_FLAGS ""
    # Use a custom .props file to set up Visual D (import paths, predefined versions...).
    #VS_USER_PROPS "${PROJECT_SOURCE_DIR}/cmake-d/VisualD.props"
)
#set_target_properties(
#  sphere_database PROPERTIES
#    RUNTIME_OUTPUT_DIRECTORY ./
#    LIBRARY_OUTPUT_DIRECTORY ./
#    ARCHIVE_OUTPUT_DIRECTORY ./
#    ARCHIVE_OUTPUT_NAME sphere_database
#    LIBRARY_OUTPUT_NAME sphere_database
#    RUNTIME_OUTPUT_NAME sphere_database
#    COMPILE_FLAGS ""
#    LINKER_LANGUAGE CXX
#    LINK_FLAGS "${SANITIZE_LDFLAGS}"
#    # Use a custom .props file to set up Visual D (import paths, predefined versions...).
#    VS_USER_PROPS "${PROJECT_SOURCE_DIR}/cmake-d/VisualD.props"
#)

add_library(interpolator interpolator.cpp)
target_include_directories(interpolator BEFORE PRIVATE
  /opt/homebrew/include/eigen3
  /usr/local/include/eigen3
  /usr/include/eigen3
)
target_link_libraries(interpolator CGAL::CGAL CGAL::CGAL_Core)

add_subdirectory(external/splinter)

#if(MSVC_IDE) # Visual Studio generator
#    # Add the .d files as (Visual D) source files to this lib, so that they show up somewhere.
#    #set(LDC_LIB_EXTRA_SOURCES ${LDC_D_SOURCE_FILES})
#    message(NOTICE "Setting D")
#    set_property(SOURCE ${CMAKE_CURRENT_LIST_DIR}/trajectory.d PROPERTY VS_TOOL_OVERRIDE "DCompile")
#    # 'Clear' the original list for the custom commands below, producing ldc2.exe and ldc2-unittest.exe -
#    # we still need a dummy .d file.
#    #set(LDC_D_SOURCE_FILES "${PROJECT_SOURCE_DIR}/dmd/root/man.d")
#    # Mark this main library target as (bold) startup project for the generated Visual Studio solution.
#    #set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${LDC_LIB})
#endif()

add_library(trajectory ${CMAKE_CURRENT_LIST_DIR}/trajectory.d wgs84tocartesian.d worldtransforms.d splines.d splinter.cpp)
set_property(TARGET trajectory PROPERTY CXX_STANDARD 17)
target_include_directories(trajectory BEFORE PRIVATE
  /opt/homebrew/include/eigen3
  /usr/local/include/eigen3
  /usr/include/eigen3
)
#Qadd_library(trajectory splinter.cpp)

set_target_properties(
    trajectory PROPERTIES
    #RUNTIME_OUTPUT_DIRECTORY ./
    #LIBRARY_OUTPUT_DIRECTORY ./
    #ARCHIVE_OUTPUT_DIRECTORY ./
    #ARCHIVE_OUTPUT_NAME trajectory
    #LIBRARY_OUTPUT_NAME trajectory
    #RUNTIME_OUTPUT_NAME trajectory
    #COMPILE_FLAGS ""
    LINKER_LANGUAGE CXX
    #LINK_FLAGS ""
    # Use a custom .props file to set up Visual D (import paths, predefined versions...).
    #VS_USER_PROPS "${PROJECT_SOURCE_DIR}/cmake-d/VisualD.props"
)

target_link_libraries(trajectory splinter-static)

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
#set(views_dir ${CMAKE_CURRENT_SOURCE_DIR}/external/matplotlibd/views)
#set(LDC_LINK_MANUALLY TRUE)

#build_d_executable(
#  niceops
#  ./niceops
#  ${CMAKE_CURRENT_SOURCE_DIR}/niceops.d
#  "-J${views_dir} -Iexternal/matplotlibd/source -Iexternal/numd/source"
#  "-Wl,--copy-dt-needed-entries"
#  ""
#  ""
#  "sphere_database;interpolator;trajectory;matplotlib-d;numd"
#  TRUE
#)
add_executable(niceops niceops.d)
target_link_libraries(niceops sphere_database interpolator trajectory matplotlib-d numd)

add_executable(wopwop2aam wopwop2aam.d)
target_link_libraries(wopwop2aam sphere_database interpolator trajectory matplotlib-d numd)

set(views_dir ${CMAKE_CURRENT_SOURCE_DIR}/external/matplotlibd/views)
target_compile_options(niceops PRIVATE "-J ${views_dir}")

#set(views_dir ${CMAKE_CURRENT_SOURCE_DIR}/external/matplotlibd/views)
target_compile_options(wopwop2aam PRIVATE "-J ${views_dir}")

if(MSVC)
  set_target_properties(
    niceops PROPERTIES
      #RUNTIME_OUTPUT_DIRECTORY ./
      #LIBRARY_OUTPUT_DIRECTORY ./
      #ARCHIVE_OUTPUT_DIRECTORY ./
      #ARCHIVE_OUTPUT_NAME sphere_database
      #LIBRARY_OUTPUT_NAME sphere_database
      #RUNTIME_OUTPUT_NAME sphere_database
      #COMPILE_FLAGS ""
      LINKER_LANGUAGE D
      LINK_FLAGS "-L/NODEFAULTLIB:MSVCRT -L/NODEFAULTLIB:libucrt ucrt.lib"
      # Use a custom .props file to set up Visual D (import paths, predefined versions...).
      #VS_USER_PROPS "${PROJECT_SOURCE_DIR}/cmake-d/VisualD.props"
  )
else()
  set_target_properties(
    niceops PROPERTIES
      #RUNTIME_OUTPUT_DIRECTORY ./
      #LIBRARY_OUTPUT_DIRECTORY ./
      #ARCHIVE_OUTPUT_DIRECTORY ./
      #ARCHIVE_OUTPUT_NAME sphere_database
      #LIBRARY_OUTPUT_NAME sphere_database
      #RUNTIME_OUTPUT_NAME sphere_database
      #COMPILE_FLAGS ""
      LINKER_LANGUAGE D
      LINK_FLAGS ""
      # Use a custom .props file to set up Visual D (import paths, predefined versions...).
      #VS_USER_PROPS "${PROJECT_SOURCE_DIR}/cmake-d/VisualD.props"
  )
endif()
add_executable(niceops_unit_tests niceops.d sphere_database.d trajectory.d splines.d)
target_link_libraries(niceops_unit_tests sphere_database interpolator trajectory matplotlib-d numd)
target_compile_options(niceops_unit_tests PRIVATE "-J ${views_dir}" )
target_compile_options(niceops_unit_tests PRIVATE "-unittest" )


# Add testing with CTest
enable_testing()
add_test(NAME regression_default COMMAND niceops  "-h" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME broken_self_test COMMAND niceops "--self_test" WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
add_test(NAME unit_test COMMAND niceops_unit_tests WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR})
