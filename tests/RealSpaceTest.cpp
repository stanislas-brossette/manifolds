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

#include <manifolds/defs.h>
#include <manifolds/utils.h>
#include <manifolds/mnf_assert.h>
#include <manifolds/RealSpace.h>
#include <manifolds/Point.h>

#ifndef _WIN32
#define BOOST_TEST_MODULE Manifolds 
#endif

#include <boost/test/unit_test.hpp>

using namespace mnf;

BOOST_AUTO_TEST_CASE(RealSpaceConstructor)
{
  RealSpace R5(5);
  BOOST_CHECK_EQUAL(R5.dim(), 5);
  BOOST_CHECK_EQUAL(R5.representationDim(), 5);
  BOOST_CHECK_EQUAL(R5.numberOfSubmanifolds(), 1);
  BOOST_CHECK(R5.isElementary());
}

BOOST_AUTO_TEST_CASE(RealPointConstructor)
{
  RealSpace R3(3);
  Point x = R3.createPoint();
  Eigen::Vector3d vy;
  vy << 1,2,3;
  Point y = R3.createPoint(vy);
  BOOST_CHECK_EQUAL(x.value().size(), 3);
  BOOST_CHECK_EQUAL(y.value().size(), 3);
  BOOST_CHECK_EQUAL(y.value()[0], 1);
  BOOST_CHECK_EQUAL(y.value()[1], 2);
  BOOST_CHECK_EQUAL(y.value()[2], 3);
}

BOOST_AUTO_TEST_CASE(RandomRealPointConstructor)
{
  RealSpace R3(3);
  Point y = R3.createRandomPoint();
  BOOST_CHECK_EQUAL(y.value().size(), 3);
  BOOST_CHECK(y.isInM());
}

BOOST_AUTO_TEST_CASE(RealSpaceZero)
{
  RealSpace R3(3);
  Point x = R3.getZero();
  for(long i = 0; i < x.value().size(); ++i)
  {
    BOOST_CHECK_EQUAL(x.value()[i], 0);
  }
}

BOOST_AUTO_TEST_CASE(RealPointIncrement)
{
  RealSpace R3(3);
  Eigen::Vector3d vx;
  Eigen::Vector3d vy;
  vx << -7, 2, 1.2;
  vy << 1,2,3;
  R3.retractation(vx,vx,vy);
  BOOST_CHECK_EQUAL(vx.size(), 3);
  BOOST_CHECK_EQUAL(vx[0], -6);
  BOOST_CHECK_EQUAL(vx[1], 4);
  BOOST_CHECK_EQUAL(vx[2], 4.2);
}

BOOST_AUTO_TEST_CASE(RealPointAddition)
{
  RealSpace R3(3);
  Eigen::Vector3d y = Eigen::Vector3d::Zero();
  Eigen::Vector3d v;
  v << 1,2,3;
  R3.retractation(y,y,v);
  R3.retractation(y,y,v);
  BOOST_CHECK_EQUAL(y.size(), 3);
  BOOST_CHECK_EQUAL(y[0], 2);
  BOOST_CHECK_EQUAL(y[1], 4);
  BOOST_CHECK_EQUAL(y[2], 6);
}

BOOST_AUTO_TEST_CASE(RealPointSubstraction)
{
  RealSpace R3(3);
  Eigen::Vector3d x;
  Eigen::Vector3d y;
  x << 4,3.4,7;
  y << 1,2,3;
  Eigen::Vector3d z; 
  R3.pseudoLog(z,x,y);
  BOOST_CHECK_EQUAL(z[0], -3);
  BOOST_CHECK_EQUAL(z[1], -1.4);
  BOOST_CHECK_EQUAL(z[2], -4);
}

BOOST_AUTO_TEST_CASE(RealPointpseudoLog0)
{
  RealSpace Space(7);
  Eigen::VectorXd x = Space.getZero().value();
  Eigen::VectorXd vy = Eigen::VectorXd::Random(Space.dim());;
  Space.retractation(x, x, vy);
  Eigen::VectorXd z(Space.dim());
  Space.pseudoLog0(z, x); 
  BOOST_CHECK(z.isApprox(vy));
}

BOOST_AUTO_TEST_CASE(RealPointDiff)
{
  RealSpace R7(7);
  Eigen::MatrixXd J;
  Eigen::VectorXd x = R7.createPoint().value();
  J = R7.diffRetractation(x);
  BOOST_CHECK(J.isIdentity());
}

BOOST_AUTO_TEST_CASE(RealApplyDiff)
{
  int c = 5;
  RealSpace Space(7);
  Index dim = Space.dim();
  Index repDim = Space.representationDim();
  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(c,repDim);
  Eigen::VectorXd x = Space.getZero().value();
  Space.retractation(x,x,Eigen::VectorXd::Random(dim));
  Eigen::MatrixXd expectedRes;
  expectedRes = Jf*Space.diffRetractation(x);
  Eigen::MatrixXd J(c,dim);
  Space.applyDiffRetractation(J, Jf, x);
  BOOST_CHECK(expectedRes.isApprox(J));
}

BOOST_AUTO_TEST_CASE(RealApplyInvDiff)
{
  int c = 5;
  RealSpace Space(7);
  Index dim = Space.dim();
  Index repDim = Space.representationDim();
  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(c,dim);
  Eigen::VectorXd x = Space.getZero().value();
  Space.retractation(x, x, Eigen::VectorXd::Random(dim));
  Eigen::MatrixXd expectedRes;
  expectedRes = Jf*Space.diffPseudoLog0(x);
  Eigen::MatrixXd J(c,repDim);
  Space.applyDiffPseudoLog0(J, Jf, x);
  BOOST_CHECK(expectedRes.isApprox(J));
}

BOOST_AUTO_TEST_CASE(RealTransport)
{
  int c = 5;
  RealSpace Space(7);
  Index dim = Space.dim();
  Eigen::MatrixXd H = Eigen::MatrixXd::Random(dim,c);
  Eigen::MatrixXd Hout(dim,c);
  Eigen::VectorXd v = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd x = Space.getZero().value();
  Space.retractation(x, x, v);
  Eigen::MatrixXd expectedRes = H;
  Space.applyTransport(Hout, H, x, v);
  BOOST_CHECK(expectedRes.isApprox(Hout));
}

BOOST_AUTO_TEST_CASE(RealInvTransport)
{
  int c = 3;
  RealSpace Space(9);
  Index dim = Space.dim();
  Eigen::MatrixXd H = Eigen::MatrixXd::Random(c,dim);
  Eigen::MatrixXd Hout(c,dim);
  Eigen::VectorXd v = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd x = Space.getZero().value();
  Space.retractation(x, x, v);
  Eigen::MatrixXd expectedRes = H;
  Space.applyInvTransport(Hout, H, x, v);
  BOOST_CHECK(expectedRes.isApprox(Hout));
}

BOOST_AUTO_TEST_CASE(RealLimitMap)
{
  RealSpace Space(9);
  Index dim = Space.dim();
  Eigen::VectorXd res(dim);
  Space.limitMap(res);
  Eigen::VectorXd expectedRes(dim);
  double i = std::numeric_limits<double>::infinity();
  expectedRes << i,i,i,i,i,i,i,i,i;
  BOOST_CHECK_EQUAL(expectedRes, res);
}

#if   EIGEN_WORLD_VERSION > 3 \
  || (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION > 2) \
  || (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION == 2 && EIGEN_MINOR_VERSION > 0)
BOOST_AUTO_TEST_CASE(RealNoAllocation)
{
  //We only test here that the operations on the manifold do not create
  //temporary. Passing arguments that are not recognize by the Eigen::Ref will
  //create temporaries, but this is the user's fault.
  RealSpace R(4);
  Index c = 3;
  Eigen::VectorXd x = Eigen::VectorXd::Random(R.representationDim());
  Eigen::VectorXd y = Eigen::VectorXd::Random(R.representationDim());
  Eigen::VectorXd v = Eigen::VectorXd::Random(R.dim());
  Eigen::VectorXd z(R.representationDim());
  Eigen::MatrixXd J0 = Eigen::MatrixXd::Random(c, R.representationDim());
  Eigen::MatrixXd J1(c, R.representationDim());
  Eigen::MatrixXd J2(c, R.representationDim());
  Eigen::MatrixXd H0 = Eigen::MatrixXd::Random(R.dim(), R.dim());
  Eigen::MatrixXd H1 = Eigen::MatrixXd::Random(R.dim(), R.dim());
  Eigen::MatrixXd H2 = Eigen::MatrixXd::Random(R.dim(), R.dim());

  Eigen::internal::set_is_malloc_allowed(false);
  utils::set_is_malloc_allowed(false);
  {
    R.retractation(z, x, v);
    R.pseudoLog(z, x, y);
    R.pseudoLog0(z, x);
    R.applyDiffRetractation(J1, J0, x);
    R.applyDiffPseudoLog0(J2, J0, x);
    R.applyTransport(H1, H0, x, v);
    R.applyInvTransport(H2, H0, x, v);
  }
  utils::set_is_malloc_allowed(true);
  Eigen::internal::set_is_malloc_allowed(true);
}
#endif
