INCLUDE_DIRECTORIES(BEFORE)

FIND_PACKAGE(Boost REQUIRED COMPONENTS unit_test_framework)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN -DPGS_ASSERT_THROW)

add_executable(RealSpaceTest RealSpaceTest.cpp)
target_link_libraries(RealSpaceTest manifoldsTest ${Boost_LIBRARIES})
add_test(RealSpaceTest RealSpaceTest)

add_executable(SO3MatrixTest SO3MatrixTest.cpp)
target_link_libraries(SO3MatrixTest manifoldsTest ${Boost_LIBRARIES})
add_test(SO3MatrixTest SO3MatrixTest)

add_executable(SO3QuaternionTest SO3QuaternionTest.cpp)
target_link_libraries(SO3QuaternionTest manifoldsTest ${Boost_LIBRARIES})
add_test(SO3QuaternionTest SO3QuaternionTest)

add_executable(S2Test S2Test.cpp)
target_link_libraries(S2Test manifoldsTest ${Boost_LIBRARIES})
add_test(S2Test S2Test)

add_executable(CartesianProductTest CartesianProductTest.cpp)
target_link_libraries(CartesianProductTest manifoldsTest ${Boost_LIBRARIES})
add_test(CartesianProductTest CartesianProductTest)

add_executable(PointTest PointTest.cpp)
target_link_libraries(PointTest manifoldsTest ${Boost_LIBRARIES})
add_test(PointTest PointTest)

add_executable(ManifoldTest ManifoldTest.cpp)
target_link_libraries(ManifoldTest manifoldsTest ${Boost_LIBRARIES})
add_test(ManifoldTest ManifoldTest)

add_executable(Manifail Manifail.cpp)
target_link_libraries(Manifail manifolds)

add_executable(TypeIdTest typeId.cpp)
target_link_libraries(TypeIdTest manifoldsTest ${Boost_LIBRARIES})
add_test(TypeIdTest TypeIdTest)