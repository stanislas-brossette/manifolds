INCLUDE_DIRECTORIES(BEFORE)

FIND_PACKAGE(Boost REQUIRED COMPONENTS unit_test_framework)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

add_executable(PGSolverExeTest test1.cpp)
target_link_libraries(PGSolverExeTest PGSolver ${Boost_LIBRARIES})

add_test(PGSolverTest PGSolverExeTest)
