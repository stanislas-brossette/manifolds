INCLUDE_DIRECTORIES(BEFORE)

FIND_PACKAGE(Boost REQUIRED COMPONENTS unit_test_framework)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

add_executable(RealSpaceTest RealSpaceTest.cpp)
target_link_libraries(RealSpaceTest PGSolverTest ${Boost_LIBRARIES})
add_test(RealSpaceTest RealSpaceTest)

add_executable(SO3Test SO3Test.cpp)
target_link_libraries(SO3Test PGSolverTest ${Boost_LIBRARIES})
add_test(SO3Test SO3Test)

add_executable(CartesianProductTest CartesianProductTest.cpp)
target_link_libraries(CartesianProductTest PGSolverTest ${Boost_LIBRARIES})
add_test(CartesianProductTest CartesianProductTest)

add_executable(testTest testTest.cpp)
target_link_libraries(testTest PGSolverTest ${Boost_LIBRARIES})
add_test(testTest testTest)
