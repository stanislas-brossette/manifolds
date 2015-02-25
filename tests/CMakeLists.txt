INCLUDE_DIRECTORIES(BEFORE)

FIND_PACKAGE(Boost REQUIRED COMPONENTS unit_test_framework)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN -DPGS_ASSERT_THROW)

add_executable(RealSpaceTest RealSpaceTest.cpp)
target_link_libraries(RealSpaceTest manifoldsTest ${Boost_LIBRARIES})
add_test(RealSpaceTest RealSpaceTest)

add_executable(SO3Test SO3Test.cpp)
target_link_libraries(SO3Test manifoldsTest ${Boost_LIBRARIES})
add_test(SO3Test SO3Test)

add_executable(CartesianProductTest CartesianProductTest.cpp)
target_link_libraries(CartesianProductTest manifoldsTest ${Boost_LIBRARIES})
add_test(CartesianProductTest CartesianProductTest)

add_executable(PointTest PointTest.cpp)
target_link_libraries(PointTest manifoldsTest ${Boost_LIBRARIES})
add_test(PointTest PointTest)

add_executable(ManifoldTest ManifoldTest.cpp)
target_link_libraries(ManifoldTest manifoldsTest ${Boost_LIBRARIES})
add_test(ManifoldTest ManifoldTest)
