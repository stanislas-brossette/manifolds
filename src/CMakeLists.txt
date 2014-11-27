set(SOURCES main.cpp)
set(HEADERS Point.hpp PointOfR.hpp)

add_library(PGSolver SHARED ${SOURCES} ${HEADERS})

set_target_properties(PGSolver PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION 0.1
)

add_executable(PGSolverExe main.cpp)
install(TARGETS PGSolver
  DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})
