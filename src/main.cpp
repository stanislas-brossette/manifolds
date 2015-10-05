// Copyright (c) 2015 CNRS
// Authors: Stanislas Brossette, Adrien Escande

// This file is part of manifolds
// manifolds is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.

// manifolds is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// manifolds. If not, see
// <http://www.gnu.org/licenses/>.

#include <iostream>
#include <stdexcept>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include <Eigen/Core>

#include <manifolds/mnf_assert.h>
#include <manifolds/RealSpace_Base.h>
#include <manifolds/SO3_Base.h>
#include <manifolds/CartesianProduct_Base.h>
#include <manifolds/Point.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/ReusableTemporaryMap.h>

using namespace mnf;

void LifeAndDeathOfAPoint()
{
  RealSpace_Base R3(3);
  Point x = R3.createPoint();
}

int main()
{
  LifeAndDeathOfAPoint();
/*srand((unsigned)time(NULL));
std::cout << "Using: Eigen" << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION
<<"." << EIGEN_MINOR_VERSION<< std::endl;
RealSpace R2(2);
RealSpace R3(3);
SO3<ExpMapMatrix> RotSpace;
CartesianProduct R2R3R2(R2, R3);
R2R3R2.multiply(R2);
CartesianProduct SO3R2R3R2(RotSpace, R2R3R2);
CartesianProduct SO3R2R3R2SO3(SO3R2R3R2, RotSpace);
CartesianProduct R3SO3(R3, RotSpace);
CartesianProduct R3SO3R3SO3(R3, RotSpace);
R3SO3R3SO3.multiply(R3);
R3SO3R3SO3.multiply(RotSpace);*/
#ifdef _WIN32
  system("pause");
#endif
  return 0;
}

