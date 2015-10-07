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
#include <manifolds/RealSpace.h>
#include <manifolds/SO3.h>
#include <manifolds/S2.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/Point.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/ReusableTemporaryMap.h>

using namespace mnf;

void LifeAndDeathOfAPoint()
{
  RealSpace R3(3);
  Point x = R3.createPoint();
}

int main()
{
srand((unsigned)time(NULL));
std::cout << "Using: Eigen" << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION
<<"." << EIGEN_MINOR_VERSION<< std::endl;
std::cout << "RealSpace R0(2);" << std::endl;
RealSpace R0(2);
std::cout << "RealSpace R1 = R0.deepCopy();" << std::endl;
RealSpace R1 = R0.deepCopy();
std::cout << "RealSpace R2(2);" << std::endl;
RealSpace R2(2);
std::cout << "R2.name()=" << R2.name() << std::endl;
std::cout << "R2.ptr()->name()=" << R2.ptr()->name() << std::endl;
std::cout << "R2.dim()=" << R2.dim() << std::endl;

std::cout << "R2.deepCopy().name() = " << R2.deepCopy().name() << std::endl;
std::cout << "R2.shallowCopy().name() = " << R2.shallowCopy().name() << std::endl;

RealSpace R3(3);
std::cout << "R3.name()=" << R3.name() << std::endl;
std::cout << "R3.ptr()->name()=" << R3.ptr()->name() << std::endl;
std::cout << "R3.dim()=" << R3.dim() << std::endl;
S2 s2;
std::cout << "s2.name()=" << s2.name() << std::endl;
std::cout << "s2.ptr()->name()=" << s2.ptr()->name() << std::endl;
std::cout << "s2.dim()=" << s2.dim() << std::endl;
SO3<ExpMapMatrix> so3;
std::cout << "so3.name()=" << so3.name() << std::endl;
std::cout << "so3.ptr()->name()=" << so3.ptr()->name() << std::endl;
std::cout << "so3.dim()=" << so3.dim() << std::endl;
CartesianProduct R2R3(R2, R3);
std::cout << "R2R3.name()=" << R2R3.name() << std::endl;
std::cout << "R2R3.ptr()->name()=" << R2R3.ptr()->name() << std::endl;
std::cout << "R2R3.dim()=" << R2R3.dim() << std::endl;
std::cout << "R2R3(0).ptr()->name()=" << R2R3(0).ptr()->name() << std::endl;
std::cout << "R2R3(0).name()=" << R2R3(0).name() << std::endl;
std::cout << "R2R3(0).dim()=" << R2R3(0).dim() << std::endl;
std::cout << "R2R3(1).ptr()->name()=" << R2R3(1).ptr()->name() << std::endl;
std::cout << "R2R3(1).name()=" << R2R3(1).name() << std::endl;
std::cout << "R2R3(1).dim()=" << R2R3(1).dim() << std::endl;
std::cout << "R2R3.deepCopy().name() = " << R2R3.deepCopy().name() << std::endl;
std::cout << "R2R3.shallowCopy().name() = " << R2R3.shallowCopy().name() << std::endl;
CartesianProduct R2R3R2(R2R3, R2);
std::cout << "R2R3R2.name()=" << R2R3R2.name() << std::endl;
std::cout << "R2R3R2.ptr()->name()=" << R2R3R2.ptr()->name() << std::endl;
std::cout << "R2R3R2.dim()=" << R2R3R2.dim() << std::endl;
//SO3<ExpMapMatrix> RotSpace;
//CartesianProduct SO3R2R3R2(RotSpace, R2R3R2);
//CartesianProduct SO3R2R3R2SO3(SO3R2R3R2, RotSpace);
//CartesianProduct R3SO3(R3, RotSpace);
//CartesianProduct R3SO3R3SO3(R3, RotSpace);
//R3SO3R3SO3.multiply(R3);
//R3SO3R3SO3.multiply(RotSpace);
//std::cout << R3SO3R3SO3.name() << std::endl;
//LifeAndDeathOfAPoint();
#ifdef _WIN32
  system("pause");
#endif
  return 0;
}

