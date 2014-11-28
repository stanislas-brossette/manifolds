set(SOURCES main.cpp CartesianProduct.cpp Manifold.cpp Point.cpp RealSpace.cpp
  SO3.cpp)
set(HEADERS CartesianProduct.h defs.h Manifold.h Point.h RealSpace.h SO3.h)

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
