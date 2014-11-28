set(SOURCES main.cpp CartesianProduct.cpp Manifold.cpp Point.cpp RealSpace.cpp
  SO3.cpp)
set(HEADERS
  ../include/pgsolver/CartesianProduct.h
  ../include/pgsolver/defs.h
  ../include/pgsolver/Manifold.h
  ../include/pgsolver/Point.h
  ../include/pgsolver/RealSpace.h
  ../include/pgsolver/SO3.h
)

add_library(PGSolver SHARED ${SOURCES} ${HEADERS})

set_target_properties(PGSolver PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION 0.1
)

target_link_libraries(PGSolver)

install(TARGETS PGSolver
  DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})
add_executable(PGSolverExe main.cpp)
target_link_libraries(PGSolverExe PGSolver)
