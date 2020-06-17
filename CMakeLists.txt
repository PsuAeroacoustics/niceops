cmake_minimum_required(VERSION 3.1)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake-d)

project(NICEOPS VERSION 1.0 LANGUAGES C CXX D)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX
    "build" CACHE PATH "default install prefix" FORCE
    )
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

#Set dylib paths relative to executable path under OSX
set(CMAKE_MACOSX_RPATH 0)
set(CMAKE_INSTALL_NAME_DIR @executable_path)
#Set attributes for Xcode project
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++14")
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
set(CMAKE_CXX_STANDARD 14)

include_directories(SYSTEM /usr/local/include /opt/local/include /sw/include /usr/include external/splinter/include external/matplotlibd/source /usr/local/include/eigen3 /usr/include/eigen3)
# Including the source directory is apparently required to support out of source builds of D modules
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_definitions(-D_THREAD_SAFE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
          "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
          FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Force shared version to be used to ensure good linking . (Some systems like Homebrew install static and shared in the same place)
find_file(NCDF NAMES libnetcdf.dylib libnetcdf.so libnetcdf.dll PATHS /usr/local/lib/ /opt/local/lib/ /sw/lib/ /usr/lib64/ /usr/lib/)

find_package(Boost)
find_package(CGAL COMPONENTS Core)
include(${CGAL_USE_FILE})
find_package (Eigen3 REQUIRED NO_MODULE)


add_library(sphere_database sphere_database.d netcdf.d)
target_link_libraries(sphere_database ${NCDF})
add_library(interpolator interpolator.cpp)
target_link_libraries(interpolator ${CGAL})

add_subdirectory(external/splinter)
add_library(trajectory trajectory.d splines.d splinter.cpp)
target_link_libraries(trajectory splinter-static)

add_subdirectory(external/matplotlibd)
add_executable(niceops niceops.d)
target_link_libraries(niceops sphere_database interpolator trajectory matplotlib-d)
set(views_dir ${CMAKE_CURRENT_SOURCE_DIR}/external/matplotlibd/views)
target_compile_options(niceops PRIVATE "-J ${views_dir}")