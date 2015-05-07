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

#ifndef _WIN32
#define BOOST_TEST_MODULE Manifold 
#endif

#include <boost/test/unit_test.hpp>

#include <manifolds/defs.h>
#include <manifolds/utils.h>
#include <manifolds/Point.h>
#include <manifolds/RealSpace.h>
#include <manifolds/SO3.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/CartesianProduct.h>

using namespace mnf;

//
//All the calculations are tested in the manifolds. Here we only test that the
//points methods are correctly mapped onto the manifolds methods
//

BOOST_AUTO_TEST_CASE(PointConstructor)
{
  RealSpace R3(3);
  Eigen::Vector3d v(1,2,3);
  const Point x = R3.createPoint(v);
  x.value();
  BOOST_CHECK(true/*v.isApprox(x.value())*/);
}

BOOST_AUTO_TEST_CASE(PointIncrement)
{
  RealSpace R3(3);
  Eigen::Vector3d v(1,2,3);
  Point x = R3.getZero();
  x.increment(v);
  BOOST_CHECK(v.isApprox(x.value()));
}

BOOST_AUTO_TEST_CASE(PointRetractation)
{
  RealSpace R3(3);
  Eigen::Vector3d v(1,2,3);
  Point x = R3.getZero();
  Point y = x.retractation(v);
  BOOST_CHECK(v.isApprox(y.value()));
}

BOOST_AUTO_TEST_CASE(PointAccessors)
{
  RealSpace R3(3);
  Eigen::Vector3d v3(1,2,3);
  RealSpace R2(2);
  Eigen::Vector2d v2(10,20);
  Eigen::VectorXd v(5);
  v << 1,2,3,10,20;
  CartesianProduct S(R3, R2);
  Point x = S.createPoint(v);
  BOOST_CHECK(v3.isApprox(x(0).value()));
  BOOST_CHECK(v2.isApprox(x(1).value()));
  BOOST_CHECK(v3.isApprox(x[0]));
  BOOST_CHECK(v2.isApprox(x[1]));
}

BOOST_AUTO_TEST_CASE(PointDimensions)
{
  SO3<ExpMapMatrix> S;
  Point p = S.createRandomPoint();
  BOOST_CHECK_EQUAL(p.getDimM(),3);
  BOOST_CHECK_EQUAL(p.getTangentDimM(),3);
  BOOST_CHECK_EQUAL(p.getRepresentationDimM(),9);
}

BOOST_AUTO_TEST_CASE(PointAddition)
{
  RealSpace R3(3);
  Eigen::Vector3d v(1,2,3);

  Point x = R3.createPoint(v);
  x = x + v;
  BOOST_CHECK(x.value().isApprox(2*v));
}


BOOST_AUTO_TEST_CASE(PointPseudoLog)
{
  RealSpace R3(3);
  Eigen::Vector3d v(1,2,3);
  Point x = R3.createPoint(v);
  Point y = R3.createPoint(2*v);
  Eigen::Vector3d res;
  res.setZero();
  res = y.pseudoLog(x);
  BOOST_CHECK(res.isApprox(-v));
}
BOOST_AUTO_TEST_CASE(PointPseudoLog0)
{
  RealSpace R3(3);
  Eigen::Vector3d v(1,2,3);

  Point y = R3.createPoint(2*v);
  Eigen::Vector3d res;

  res = y.pseudoLog0();
  BOOST_CHECK(res.isApprox(2*v));
}
BOOST_AUTO_TEST_CASE(PointMinus)
{
  RealSpace R3(3);
  Eigen::Vector3d v(1,2,3);

  Point x = R3.createPoint(v);
  Point y = R3.createPoint(2*v);
  Eigen::Vector3d res;

  res = y - x;
  BOOST_CHECK(res.isApprox(v));
}

BOOST_AUTO_TEST_CASE(SubPointManipulation)
{
  RealSpace R2(2);
  RealSpace R3(3);
  RealSpace R5(5);
  RealSpace R8(8);
  RealSpace R13(13);

  /*
              P
            /   \
           /      P3
          /      /  \
        P1      P2  R13
       /  \    /  \
      R2  R3  P5  P8
  */
  CartesianProduct P1(R2, R3);
  CartesianProduct P2(R5, R8);
  CartesianProduct P3(P2, R13);
  CartesianProduct P(P1, P3);
  Eigen::VectorXd v = Eigen::VectorXd::LinSpaced(P.representationDim(), 1, static_cast<double>(P.representationDim()));
  Point x = P.createPoint(v);

  ConstSubPoint p1 = x(0);  //const subpoint
  SubPoint p3 = x(1);       //non const subpoint
  SubPoint p2 = p3(0);      //non const again
  SubPoint r2 = x(0)(0);    //chaining subpoint
  ConstSubPoint r3 = p1(1); //const subpoint of a const subpoint
  //SubPoint r3b = p1(1);   //non const subpoint to a const subpoint. Does not compile, normal
  SubPoint r5 = p2(0);      //non const subpoint to a non const subpoint
  SubPoint r8 = p2(1);      //idem
  ConstSubPoint r13 = p3(1);//const subpoint of a non-const subpoint

  //Correct submanifold
  BOOST_CHECK(&x.getManifold() == &P);
  BOOST_CHECK(&r2.getManifold() == &R2);
  BOOST_CHECK(&r3.getManifold() == &R3);
  BOOST_CHECK(&r5.getManifold() == &R5);
  BOOST_CHECK(&r8.getManifold() == &R8);
  BOOST_CHECK(&r13.getManifold() == &R13);
  BOOST_CHECK(&p1.getManifold() == &P1);
  BOOST_CHECK(&p2.getManifold() == &P2);
  BOOST_CHECK(&p3.getManifold() == &P3);

  //Changing values of a subpoint
  Eigen::VectorXd v5(5); v5 << 0, 0, 7, 0, 0;
  r5.value() = v5;
  BOOST_CHECK(p2.value().head<5>() == v5);
  BOOST_CHECK(x.value().segment<5>(5) == v5);

  //New point is a duplicate, not a subpoint
  Eigen::VectorXd w5(5); w5 << 7, 7, 0, 7, 7;
  Point newr5 = p2(0);
  newr5.value() = w5;
  BOOST_CHECK(p2.value().head<5>() == v5);
  BOOST_CHECK(x.value().segment<5>(5) == v5);
}

BOOST_AUTO_TEST_CASE(PointDiffRetractation)
{
  SO3<ExpMapMatrix> S;
  Eigen::MatrixXd J;
  Eigen::MatrixXd Jtest(9,3);
  Jtest <<         0,  0.003580526006716, -0.558259982176135,
                   0,  0.646068748944272,  0.634553549147103,
                   0, -0.763270824459509,  0.534497507539106,
  -0.003580526006716,                  0, -0.829658346630838,
  -0.646068748944272,                  0, -0.424189774632061,
   0.763270824459509,                  0, -0.362946363755562,
   0.558259982176135,  0.829658346630838,                  0,
  -0.634553549147103,  0.424189774632061,                  0,
  -0.534497507539106,  0.362946363755562,                  0;
  Eigen::Vector3d v(0.680375, -0.211234, 0.566198);
  Point x = S.getZero();
  x.retractation(x, v);
  J = x.diffRetractation();
  BOOST_CHECK(J.isApprox(Jtest));
}

BOOST_AUTO_TEST_CASE(PointApplyDiff)
{
  int c = 5;
  SO3<ExpMapMatrix> S;
  Index dim = S.dim();
  Index repDim = S.representationDim();
  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(c,repDim);
  Point x = S.getZero();
  x.retractation(x, Eigen::VectorXd::Random(dim));
  Eigen::MatrixXd expectedRes;
  expectedRes = Jf*x.diffRetractation();
  Eigen::MatrixXd J(c,dim);
  x.applyDiffRetractation(J, Jf);
  BOOST_CHECK(expectedRes.isApprox(J));
}

BOOST_AUTO_TEST_CASE(PointInvDiff)
{
  SO3<ExpMapMatrix> S;
  Eigen::MatrixXd J;
  Eigen::MatrixXd Jtest(3,9);
  Jtest <<
  -0.064043491813865,                 0,                  0,                  0, -0.064043491813865, 0.545030346992499,                 0, -0.545030346992499, -0.064043491813865,
  -0.110117993664377,                 0, -0.545030346992499,                  0, -0.110117993664377,                 0, 0.545030346992499,                  0, -0.110117993664377,
  -0.042109988599266, 0.545030346992499,                  0, -0.545030346992499, -0.042109988599266,                 0,                 0,                  0, -0.042109988599266;
  Eigen::Vector3d v(0.3403857, 0.58526775, 0.223811);
  Point x = S.getZero();
  x.retractation(x,v);
  J = x.diffPseudoLog0();
  BOOST_CHECK(J.isApprox(Jtest));
}

BOOST_AUTO_TEST_CASE(SO3ApplyInvDiff)
{
  int c = 5;
  SO3<ExpMapMatrix> S;
  Index dim = S.dim();
  Index repDim = S.representationDim();
  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(c,dim);
  Point x = S.getZero();
  x.retractation(x, Eigen::VectorXd::Random(dim));
  Eigen::MatrixXd expectedRes;
  expectedRes = Jf*x.diffPseudoLog0();
  Eigen::MatrixXd J(c,repDim);
  x.applyDiffPseudoLog0(J, Jf);
  BOOST_CHECK(expectedRes.isApprox(J));
}

#if   EIGEN_WORLD_VERSION > 3 \
  || (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION > 2) \
  || (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION == 2 && EIGEN_MINOR_VERSION > 0)
BOOST_AUTO_TEST_CASE(PointNoAllocation)
{
  //We only test here that the operations on the manifold do not create
  //temporary. Passing arguments that are not recognize by the Eigen::Ref will
  //create temporaries, but this is the user's fault.
  RealSpace S(3);
  //Index dim = S.dim();
  //Index repDim = S.representationDim();

  Point x = S.createPoint();
  x.value() << 1, 2, 3;
  Point y = S.createPoint();
  y.value() << 4, 5, 6;
  Eigen::Vector3d v(9, 8, 7);

  Eigen::internal::set_is_malloc_allowed(false);
  utils::set_is_malloc_allowed(false);
  {
    std::cout << "Memory allocation tests:" << std::endl;
    x.retractation(y.value(),v);
    x.retractation(y,v);
    std::cout << "- method 'retractation' passed" << std::endl;
    x.pseudoLog(v, y);
    std::cout << "- method 'pseudoLog' passed" << std::endl;
    x.pseudoLog0(v);
    std::cout << "- method 'pseudoLog0' passed" << std::endl;
  }
  utils::set_is_malloc_allowed(true);
  Eigen::internal::set_is_malloc_allowed(true);
}

#endif

