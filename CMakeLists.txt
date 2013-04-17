CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
include(ExternalProject)

PROJECT(ssa)

# The library prefix
SET(LIB_PREFIX ssa_)

SET(ssa_C_FLAGS)
SET(ssa_CXX_FLAGS)

# default built type
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

# postfix, based on type
SET(CMAKE_DEBUG_POSTFIX "_d" CACHE STRING "postfix applied to debug build of libraries")
SET(CMAKE_RELEASE_POSTFIX "" CACHE STRING "postfix applied to release build of libraries")
SET(CMAKE_RELWITHDEBINFO_POSTFIX "_rd" CACHE STRING "postfix applied to release-with-debug-information libraries")
SET(CMAKE_MINSIZEREL_POSTFIX "_s" CACHE STRING "postfix applied to minimium-size-build libraries")

# On the Mac platform, configure the RPATH as per the install, to
# avoid the problem of loading both the built and installed versions
# of the shared targets
IF(APPLE)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
  SET(CMAKE_INSTALL_RPATH "")
ENDIF(APPLE)

# Set the output directory for the build executables and libraries
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${ssa_SOURCE_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${ssa_SOURCE_DIR}/bin)
IF(WIN32)
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ssa_SOURCE_DIR}/bin)
ELSE(WIN32)
  SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${ssa_SOURCE_DIR}/lib)
ENDIF(WIN32)

# Set search directory for looking for our custom CMake scripts to
# look for SuiteSparse and QGLViewer
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${ssa_SOURCE_DIR}/cmake_modules)

# Detect OS and define macros appropriately
IF(WIN32)
  IF(CYGWIN)
    MESSAGE(STATUS "Compiling on Cygwin")
    ADD_DEFINITIONS(-DCYGWIN)
  ELSE(CYGWIN)
    ADD_DEFINITIONS(-DWINDOWS)
    MESSAGE(STATUS "Compiling on Windows")
  ENDIF(CYGWIN)
ELSEIF(UNIX)
  ADD_DEFINITIONS(-DUNIX)
  MESSAGE(STATUS "Compiling on Unix")
ENDIF(WIN32)

SET(SSA_BUILD_SUITESPARSE OFF CACHE BOOL "Build GotoBLAS2 and Suitesparse (experimental)")
# For building the CHOLMOD / CSPARSE solvers
IF(SSA_BUILD_SUITESPARSE)
  IF(UNIX)
    SET(AMD_LIBRARY ${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/lib/libamd.a CACHE STRING "Path to amd library" FORCE)
    SET(CAMD_LIBRARY "${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/CAMD/Lib/libcamd.a" CACHE STRING "Path of camd library" FORCE)
    SET(COLAMD_LIBRARY "${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/lib/libcolamd.a" CACHE STRING "Path of colamd library" FORCE)
    SET(CCOLAMD_LIBRARY "${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/CCOLAMD/Lib/libccolamd.a" CACHE STRING "Path of ccolamd library" FORCE)

    SET(CHOLMOD_INCLUDE_DIR "${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/include" CACHE STRING "Suitesparse include dirs" FORCE)
    SET(CHOLMOD_LIBRARY "${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/CHOLMOD/Lib/libcholmod.a" CACHE STRING "Path of cholmod library" FORCE)
    SET(CSPARSE_INCLUDE_DIR "${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/include" CACHE STRING "Cxsparse include dir" FORCE) 
    SET(CSPARSE_LIBRARY "${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/lib/libcxsparse.a" CACHE STRING "Path of cxsparse library" FORCE)
    SET(BLAS_LIBRARIES "${ssa_SOURCE_DIR}/EXTERNAL/OpenBLAS2/libopenblas.so" CACHE STRING "Path of BLAS library"  FORCE)
    SET(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARY} ${AMD_LIBRARY} ${CAMD_LIBRARY} ${COLAMD_LIBRARY} ${CCOLAMD_LIBRARY} ${BLAS_LIBRARIES} CACHE STRING "Paths to cholmod" FORCE)
  ENDIF(UNIX)
  IF(WIN32) 
    SET(AMD_LIBRARY ${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/lib/libamd.lib CACHE STRING "Path to amd library" FORCE)
    SET(CAMD_LIBRARY "${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/CAMD/Lib/libcamd.lib" CACHE STRING "Path of camd library" FORCE)
    SET(COLAMD_LIBRARY "${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/lib/libcolamd.lib" CACHE STRING "Path of colamd library" FORCE)
    SET(CCOLAMD_LIBRARY "${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/CCOLAMD/Lib/libccolamd.lib" CACHE STRING "Path of ccolamd library" FORCE)

    SET(CHOLMOD_INCLUDE_DIR "${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/include" CACHE STRING "Suitesparse include dirs" FORCE)
    SET(CHOLMOD_LIBRARY "${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/CHOLMOD/Lib/libcholmod.lib" CACHE STRING "Path of cholmod library" FORCE)
    SET(CSPARSE_INCLUDE_DIR "${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/include" CACHE STRING "Cxsparse include dir" FORCE) 
    SET(CSPARSE_LIBRARY "${ssa_SOURCE_DIR}/EXTERNAL/SuiteSparse/lib/libcxsparse.lib" CACHE STRING "Path of cxsparse library" FORCE)
    SET(BLAS_LIBRARIES "${ssa_SOURCE_DIR}/EXTERNAL/OpenBLAS2/libopenblas.dll" CACHE STRING "Path of BLAS library"  FORCE)
    SET(CHOLMOD_LIBRARIES ${CHOLMOD_LIBRARY} ${AMD_LIBRARY} ${CAMD_LIBRARY} ${COLAMD_LIBRARY} ${CCOLAMD_LIBRARY} ${BLAS_LIBRARIES} CACHE STRING "Paths to cholmod" FORCE)
  ENDIF(WIN32) 

ELSE(SSA_BUILD_SUITESPARSE)
  FIND_PACKAGE(SuiteSparse)
ENDIF(SSA_BUILD_SUITESPARSE)

# Eigen library parallelise itself, though, presumably due to performance issues
# OPENMP is experimental. We experienced some slowdown with it
FIND_PACKAGE(OpenMP)
SET( SSA_USE_OPENMP OFF CACHE BOOL "Build ssa with OpenMP support (EXPERIMENTAL)")
IF(OPENMP_FOUND AND SSA_USE_OPENMP)
  SET(ssa_C_FLAGS "${ssa_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(ssa_CXX_FLAGS "${ssa_CXX_FLAGS}  ${OpenMP_CXX_FLAGS}") #-DEIGEN_DONT_PARALLELIZE
  MESSAGE(STATUS "Compiling with OpenMP support")
ENDIF(OPENMP_FOUND AND SSA_USE_OPENMP)

# FIND OR CHECKOUT G2O 
SET( SSA_USE_G2O_GIT ON CACHE BOOL "Download and build G2O from git (EXPERIMENTAL)")
IF(SSA_USE_G2O_GIT)
  MESSAGE(STATUS "Will checkout g2o git repository")
  SET(G2O_PATH ${ssa_SOURCE_DIR}/EXTERNAL/g2o CACHE PATH "Directory of G2O")
  # Set up the top-level include directories
  SET(G2O_INCLUDE_DIR ${G2O_PATH}/include CACHE PATH "Include directory of G2O")
  include_directories(${CHOLMOD_INCLUDE_DIR})
ELSE(SSA_USE_G2O_GIT)
  MESSAGE(STATUS "Search local g2o installation")
  FIND_PACKAGE(G2O)
  IF(G2O_FOUND)
    MESSAGE(STATUS "Found g2o " ${G2O_INCLUDE_DIR})
    SET(G2O_PATH ${G2O_INCLUDE_DIR}/..)
  ELSE(G2O_FOUND)
    MESSAGE(FATAL_ERROR "g2o not found! Either set SSA_USE_G2O_GIT to ON or install g2o! ")
  ENDIF(G2O_FOUND)
ENDIF(SSA_USE_G2O_GIT)

LINK_DIRECTORIES( "${G2O_PATH}/lib" )
include_directories(${G2O_INCLUDE_DIR})

FIND_PACKAGE(CUDA)
SET( SSA_USE_CUDA OFF CACHE BOOL "Build ssa with CUDA support (EXPERIMENTAL)")
IF(CUDA_FOUND AND SSA_USE_CUDA)
  SET(CUDA_SDK_LIBRARIES ${CUDA_SDK_ROOT_DIR}/shared/lib/libshrutil_x86_64.a    ${CUDA_SDK_ROOT_DIR}/C/lib/libcutil_x86_64.a  ${CUDA_SDK_ROOT_DIR}/C/lib/libparamgl_x86_64.a  ${CUDA_SDK_ROOT_DIR}/C/lib/librendercheckgl_x86_64.a)
  MESSAGE(STATUS "Compiling with CUDA support")
ENDIF(CUDA_FOUND AND SSA_USE_CUDA)

SET(SSA_USE_CUSP OFF CACHE BOOL "Build ssa with CUSP support (EXPERIMENTAL)")
IF(CUDA_FOUND AND SSA_USE_CUSP)
  SET(CUSP_INCLUDE_DIR  ${ssa_SOURCE_DIR}/EXTERNAL/cusp/ ${ssa_SOURCE_DIR}/EXTERNAL/thrust/)
ENDIF(CUDA_FOUND AND SSA_USE_CUSP)
# OpenGL is used in the draw actions for the different types, as well
# as for creating the GUI itself
# FIND_PACKAGE(OpenGL)
# SET(G2O_USE_OPENGL ON CACHE BOOL "Build g2o with OpenGL support for visualization")
# IF (OPENGL_FOUND AND G2O_USE_OPENGL)
#   SET (G2O_HAVE_OPENGL 1)
#   INCLUDE_DIRECTORIES(${OPENGL_INCLUDE})
#   MESSAGE(STATUS "Compiling with OpenGL support")
# ENDIF(OPENGL_FOUND AND G2O_USE_OPENGL)

# For building the GUI
# OpenGL is used in the draw actions for the different types, as well
# as for creating the GUI itself
FIND_PACKAGE(OpenGL)
FIND_PACKAGE(QGLViewer)

# shall we build the core apps using the library
SET(SSA_BUILD_APPS ON CACHE BOOL "Build ssa apps")
IF(SSA_BUILD_APPS)
  MESSAGE(STATUS "Compiling ssa apps")
ENDIF(SSA_BUILD_APPS)

# # shall we build the examples
# SET(G2O_BUILD_EXAMPLES ON CACHE BOOL "Build g2o examples")
# IF(G2O_BUILD_EXAMPLES)
#   MESSAGE(STATUS "Compiling g2o examples")
# ENDIF(G2O_BUILD_EXAMPLES)


find_package(PCL 1.6 REQUIRED)
include_directories(SYSTEM ${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

# Flann
#find_package(flann 1.6 REQUIRED)

# OpenNI
include_directories("/usr/include/openni")

# Compiler specific options for gcc
IF(CMAKE_COMPILER_IS_GNUCXX)
  # Generic settings for optimisation
  SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -msse2") 
  SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -msse2") 

  # OS X
  #IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    #SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}") 
    #SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}") 
  #ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # Linux
  IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native") 
    SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native") 
  ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # activate warnings !!!
  SET(ssa_C_FLAGS "${ssa_C_FLAGS} -Wall -W")
  SET(ssa_CXX_FLAGS "${ssa_CXX_FLAGS} -Wall -W")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

# specifying compiler flags
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ssa_CXX_FLAGS}")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ssa_C_FLAGS}")

include_directories(${PROJECT_BINARY_DIR} ${ssa_SOURCE_DIR})

# Include the subdirectories
ADD_SUBDIRECTORY(EXTERNAL)
ADD_SUBDIRECTORY(ssa)
ADD_SUBDIRECTORY(pcl_binding)
ADD_SUBDIRECTORY(import)
