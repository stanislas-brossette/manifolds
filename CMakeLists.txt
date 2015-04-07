cmake_minimum_required(VERSION 2.11)

if(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++0x -pedantic")
endif()

enable_testing()

set(CUSTOM_HEADER_DIR manifolds)
include(cmake/base.cmake)
include(cmake/boost.cmake)
include(cmake/eigen.cmake)
include(cmake/pkg-config.cmake)

set(PROJECT_NAME manifolds)
set(PROJECT_DESCRIPTION manifolds)
set(PROJECT_URL https://brossette@idh.lirmm.fr/././././pgsolver.git)

set(HEADERS )

setup_project()
if(CMAKE_BUILD_TYPE MATCHES Rel)
  if(NOT WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
  endif()
endif()

# Require Eigen >= 3.2.1 for Eigen::Ref
IF(NOT WIN32)
  SET(Eigen_REQUIRED "eigen3 >= 3.2.1")
  SEARCH_FOR_EIGEN()
ELSE()
  find_path(Eigen_INCLUDE_DIR NAMES signature_of_eigen3_matrix_library
    PATHS
    ${CMAKE_INSTALL_PREFIX}/include
    PATH_SUFFIXES eigen3 eigen
  )
  INCLUDE_DIRECTORIES(SYSTEM ${Eigen_INCLUDE_DIR})
ENDIF()
#IF(PKG_CONFIG_FOUND)
#ADD_REQUIRED_DEPENDENCY(EigenQP)
#ELSE()
#find_path(EigenQP_INCLUDE_DIR NAMES EigenQP/EigenQPAPI.h
#    PATHS
#    ${CMAKE_INSTALL_PREFIX}/include
#  )
#  INCLUDE_DIRECTORIES(SYSTEM ${EigenQP_INCLUDE_DIR})
#find_path(EigenQP_LIB_DIR NAMES EigenQP.lib libEigenQP.so
#    PATHS
#    ${CMAKE_INSTALL_PREFIX}/library
#    )
#   LINK_DIRECTORIES(${EigenQP_LIB_DIR})
#ENDIF()
include_directories(include)

FIND_PACKAGE(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

set(CMAKE_DEBUG_POSTFIX "_d")

# add_subdirectory(bin)
# add_subdirectory(share)
add_subdirectory(src)
add_subdirectory(tests)
# add_subdirectory(var)

setup_project_finalize()
