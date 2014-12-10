cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++0x -pedantic")

enable_testing()

include(cmake/base.cmake)
include(cmake/boost.cmake)
include(cmake/eigen.cmake)
include(cmake/pkg-config.cmake)

set(PROJECT_NAME PGSolver)
set(PROJECT_DESCRIPTION PGSolver)
set(PROJECT_URL https://brossette@idh.lirmm.fr/././././pgsolver.git)

set(CUSTOM_HEADER_DIR PGSolver)
set(HEADERS )

setup_project()
if(CMAKE_BUILD_TYPE MATCHES Rel)
  if(NOT WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
  endif()
endif()

# Add main library to pkg-config file.
PKG_CONFIG_APPEND_LIBS(PGSolver)
#
#search_for_boost()

# Require Eigen >= 3.2.1 for Eigen::Ref
SET(Eigen_REQUIRED "eigen3 >= 3.2.1")
search_for_eigen()

include_directories(include)

# add_subdirectory(bin)
# add_subdirectory(share)
add_subdirectory(src)
add_subdirectory(tests)
# add_subdirectory(var)

setup_project_finalize()
