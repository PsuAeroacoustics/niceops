cmake_minimum_required(VERSION 3.1)

project(interpolator_cpp VERSION 1.0 LANGUAGES C CXX)

message(NOTICE "BUILD_TYPE set to ${BUILD_TYPE}")
if(${BUILD_TYPE} STREQUAL "")
	set(CMAKE_BUILD_TYPE "Debug")
elseif(${BUILD_TYPE} MATCHES "debug")
	set(CMAKE_BUILD_TYPE "Debug")
elseif(${BUILD_TYPE} MATCHES "release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX
    "build" CACHE PATH "default install prefix" FORCE
    )
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

if("${CMAKE_CL_64}" STREQUAL "1")
  message(NOTICE "64 bit!")
  set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS}-m64")
endif()


#Set dylib paths relative to executable path under OSX
set(CMAKE_MACOSX_RPATH 0)
set(CMAKE_INSTALL_NAME_DIR @executable_path)
#Set attributes for Xcode project
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++17")
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
set(CMAKE_CXX_STANDARD 14)

include_directories(SYSTEM /usr/local/include /opt/local/include /sw/include /usr/include)

# Including the source directory is apparently required to support out of source builds of D modules
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
add_definitions(-D_THREAD_SAFE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
          "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
          FORCE)
endif (NOT CMAKE_BUILD_TYPE)

find_package(Boost)
find_package(CGAL COMPONENTS Core)
#include(${CGAL_USE_FILE})

if(NOT ${GMP_INCLUDES} STREQUAL "")
  include_directories(${GMP_INCLUDES})
endif()

if(NOT ${MPFR_INCLUDES} STREQUAL "")
  include_directories(${MPFR_INCLUDES})
endif()

add_library(interpolator interpolator.cpp)

set_property(TARGET interpolator PROPERTY CXX_STANDARD 20)

#target_link_libraries(interpolator ${CGAL})
target_link_libraries(interpolator CGAL::CGAL)
target_compile_options(interpolator PRIVATE -fPIE)