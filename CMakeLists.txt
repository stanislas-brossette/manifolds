cmake_minimum_required(VERSION 2.8)

#project(PGSolver)
#add_executable(PGSolver src/main.cpp)

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

# Add main library to pkg-config file.
PKG_CONFIG_APPEND_LIBS(PGSolver)
#
search_for_boost()
search_for_eigen()

# add_subdirectory(bin)
# add_subdirectory(share)
add_subdirectory(src)
# add_subdirectory(tests)
# add_subdirectory(var)

setup_project_finalize()
