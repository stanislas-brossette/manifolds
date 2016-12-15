## Copyright (c) 2015 CNRS
## Authors: Stanislas Brossette, Adrien Escande, Pierre Gergondet 
##
## This file is part of manifolds
## manifolds is free software: you can redistribute it
## and/or modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation, either version
## 3 of the License, or (at your option) any later version.
##
## manifolds is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the implied warranty
## of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Lesser Public License for more details.  You should have
## received a copy of the GNU Lesser General Public License along with
## manifolds. If not, see
## <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.8)

if(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++0x -pedantic")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++0x -pedantic -O1 -fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls")
endif()

enable_testing()

set(CUSTOM_HEADER_DIR manifolds)
include(cmake/base.cmake)
include(cmake/boost.cmake)
include(cmake/eigen.cmake)
include(cmake/pkg-config.cmake)
include(cmake/msvc-specific.cmake)

set(PROJECT_NAME manifolds)
set(PROJECT_DESCRIPTION manifolds)
set(PROJECT_URL https://brossette@idh.lirmm.fr/././././pgsolver.git)
set(PROJECT_DEBUG_POSTFIX "_d")

option(DISABLE_TESTS "Disable unit tests" OFF)

set(HEADERS )

setup_project()
if(CMAKE_BUILD_TYPE MATCHES Rel)
  if(NOT WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
  endif()
endif()

# Require Eigen >= 3.2.1 for Eigen::Ref
IF(PKG_CONFIG_FOUND)
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

# For MSVC, set local environment variable to enable finding the built dlls
# of roboptim-core and the plugins when launching ctest with RUN_TESTS
IF(MSVC)
  SET(CMAKE_MSVCIDE_RUN_PATH "\$(SolutionDir)/src/\$(Configuration)")
ENDIF(MSVC)

# add_subdirectory(bin)
# add_subdirectory(share)
add_subdirectory(src)

if(NOT ${DISABLE_TESTS})
  add_subdirectory(tests)
endif()
# add_subdirectory(var)

setup_project_finalize()
