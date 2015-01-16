set(SOURCES 
  main.cpp 
  manifolds/CartesianProduct.cpp 
  manifolds/ExpMapMatrix.cpp 
  manifolds/Manifold.cpp 
  manifolds/Point.cpp 
  manifolds/RealSpace.cpp
  manifolds/ReusableTemporaryMap.cpp
  solver/ConstraintManager.cpp
  solver/ExampleGeometricProblem.cpp
  #solver/ExampleProblem.cpp
  solver/Filter.cpp
  solver/HessianUpdater.cpp
  solver/LineSearcher.cpp
  solver/Problem.cpp
  solver/Solver.cpp
  )
set(HEADERS
  ../include/pgsolver/manifolds/CartesianProduct.h
  ../include/pgsolver/manifolds/defs.h
  ../include/pgsolver/manifolds/ExpMapMatrix.h
  ../include/pgsolver/manifolds/Manifold.h
  ../include/pgsolver/manifolds/pgs_assert.h
  ../include/pgsolver/manifolds/Point.h
  ../include/pgsolver/manifolds/RealSpace.h
  ../include/pgsolver/manifolds/ReusableTemporaryMap.h
  ../include/pgsolver/manifolds/SO3.h
  ../include/pgsolver/solver/ConstraintManager.h
  ../include/pgsolver/solver/ExampleGeometricProblem.h
  #../include/pgsolver/solver/ExampleProblem.h
  ../include/pgsolver/solver/Filter.h
  ../include/pgsolver/solver/HessianUpdater.h
  ../include/pgsolver/solver/LineSearcher.h
  ../include/pgsolver/solver/Problem.h
  ../include/pgsolver/solver/Solver.h
  ../include/pgsolver/solver/SolverOptions.h
  ../include/pgsolver/solver/TempData.h
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
  IF(PKG_CONFIG_FOUND)
  PKG_CONFIG_USE_DEPENDENCY(${lib} EigenQP)
  ELSE()
  target_link_libraries(${lib} EigenQP)
  ENDIF()
  target_compile_definitions(${lib} PUBLIC "-DPGS_EXPORT")
endforeach()

target_compile_definitions(PGSolverTest PUBLIC "-DPGS_ASSERT_THROW -DEIGEN_RUNTIME_NO_MALLOC")

install(TARGETS PGSolver
  RUNTIME DESTINATION ${CMAKE_INSTALL_FULL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}
)
if(WIN32)
install(TARGETS PGSolverTest
  RUNTIME DESTINATION ${CMAKE_INSTALL_FULL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR}
)
endif()

add_executable(PGSolverExe main.cpp)
target_link_libraries(PGSolverExe PGSolver)
