set(SOURCES 
  main.cpp 
  CartesianProduct.cpp 
  ExpMapMatrix.cpp 
  Manifold.cpp 
  Point.cpp 
  RealSpace.cpp
  ReusableTemporaryMap.cpp
  )
set(HEADERS
  ../include/pgsolver/CartesianProduct.h
  ../include/pgsolver/defs.h
  ../include/pgsolver/ExpMapMatrix.h
  ../include/pgsolver/Manifold.h
  ../include/pgsolver/pgs_assert.h
  ../include/pgsolver/Point.h
  ../include/pgsolver/RealSpace.h
  ../include/pgsolver/ReusableTemporaryMap.h
  ../include/pgsolver/SO3.h
)

set(LIBRARIES
  PGSolver
  PGSolverTest
)
foreach(lib ${LIBRARIES})
  add_library(${lib} SHARED ${SOURCES} ${HEADERS})
  set_target_properties(${lib} PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION 0.1
  )
  target_link_libraries(${lib})
endforeach()

target_compile_definitions(PGSolverTest PRIVATE "-DPGS_ASSERT_THROW -DEIGEN_RUNTIME_NO_MALLOC")

install(TARGETS PGSolver
  DESTINATION ${CMAKE_INSTALL_FULL_BINDIR})

add_executable(PGSolverExe main.cpp)
target_link_libraries(PGSolverExe PGSolver)
