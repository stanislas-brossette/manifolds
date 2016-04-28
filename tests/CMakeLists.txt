## Copyright (c) 2015 CNRS
## Authors: Stanislas Brossette, Adrien Escande 

## This file is part of manifolds
## manifolds is free software: you can redistribute it
## and/or modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation, either version
## 3 of the License, or (at your option) any later version.

## manifolds is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the implied warranty
## of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Lesser Public License for more details.  You should have
## received a copy of the GNU Lesser General Public License along with
## manifolds. If not, see
## <http://www.gnu.org/licenses/>.

INCLUDE_DIRECTORIES(BEFORE)

FIND_PACKAGE(Boost REQUIRED COMPONENTS unit_test_framework)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN -DMNF_ASSERT_THROW -DEIGEN_RUNTIME_NO_MALLOC)

MACRO(MANIFOLD_TEST NAME)
  add_executable(${NAME} ${NAME}.cpp ${ARGN})
  target_link_libraries(${NAME} manifoldsTest ${Boost_LIBRARIES})
  add_test(${NAME} ${NAME})
  # Adding a project configuration file (for MSVC only)
  GENERATE_MSVC_DOT_USER_FILE(${NAME})
ENDMACRO()

MANIFOLD_TEST(RealSpaceTest)
MANIFOLD_TEST(SO3MatrixTest)
MANIFOLD_TEST(SO3QuaternionTest)
MANIFOLD_TEST(S2Test)
MANIFOLD_TEST(CartesianProductTest)
MANIFOLD_TEST(PointTest)
MANIFOLD_TEST(ManifoldTest)
MANIFOLD_TEST(TypeIdTest)
