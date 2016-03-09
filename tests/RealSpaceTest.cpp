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
#include <manifolds/S2.h>
#include <manifolds/SO3.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/Point.h>

#ifndef _WIN32
#define BOOST_TEST_MODULE Manifolds
#endif

#include <boost/test/unit_test.hpp>

using namespace mnf;
using namespace Eigen;

BOOST_AUTO_TEST_CASE(RealSpaceConstructor)
{
  RealSpace R5(5);
  BOOST_CHECK_EQUAL(R5.dim(), 5);
  BOOST_CHECK_EQUAL(R5.representationDim(), 5);
  BOOST_CHECK_EQUAL(R5.numberOfSubManifolds(), 1);
  BOOST_CHECK(R5.isElementary());
}

BOOST_AUTO_TEST_CASE(RealPointConstructor)
{
  RealSpace R3(3);
  Point x = R3.createPoint();
  Vector3d vy;
  vy << 1, 2, 3;
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
  for (long i = 0; i < x.value().size(); ++i)
  {
    BOOST_CHECK_EQUAL(x.value()[i], 0);
  }
}

BOOST_AUTO_TEST_CASE(RealPointIncrement)
{
  RealSpace R3(3);
  Vector3d vx;
  Vector3d vy;
  vx << -7, 2, 1.2;
  vy << 1, 2, 3;
  R3.retractation(vx, vx, vy);
  BOOST_CHECK_EQUAL(vx.size(), 3);
  BOOST_CHECK_EQUAL(vx[0], -6);
  BOOST_CHECK_EQUAL(vx[1], 4);
  BOOST_CHECK_EQUAL(vx[2], 4.2);
}

BOOST_AUTO_TEST_CASE(RealPointAddition)
{
  RealSpace R3(3);
  Vector3d y = Vector3d::Zero();
  Vector3d v;
  v << 1, 2, 3;
  R3.retractation(y, y, v);
  R3.retractation(y, y, v);
  BOOST_CHECK_EQUAL(y.size(), 3);
  BOOST_CHECK_EQUAL(y[0], 2);
  BOOST_CHECK_EQUAL(y[1], 4);
  BOOST_CHECK_EQUAL(y[2], 6);
}

BOOST_AUTO_TEST_CASE(RealPointSubstraction)
{
  RealSpace R3(3);
  Vector3d x;
  Vector3d y;
  x << 4, 3.4, 7;
  y << 1, 2, 3;
  Vector3d z;
  R3.pseudoLog(z, x, y);
  BOOST_CHECK_EQUAL(z[0], -3);
  BOOST_CHECK_EQUAL(z[1], -1.4);
  BOOST_CHECK_EQUAL(z[2], -4);
}

BOOST_AUTO_TEST_CASE(RealPointpseudoLog0)
{
  RealSpace Space(7);
  VectorXd x = Space.getZero().value();
  VectorXd vy = VectorXd::Random(Space.dim());
  ;
  Space.retractation(x, x, vy);
  VectorXd z(Space.dim());
  Space.pseudoLog0(z, x);
  BOOST_CHECK(z.isApprox(vy));
}

BOOST_AUTO_TEST_CASE(RealPointDiff)
{
  RealSpace R7(7);
  MatrixXd J;
  VectorXd x = R7.createPoint().value();
  J = R7.diffRetractation(x);
  BOOST_CHECK(J.isIdentity());
}

BOOST_AUTO_TEST_CASE(RealApplyDiff)
{
  int c = 5;
  RealSpace Space(7);
  Index dim = Space.dim();
  Index repDim = Space.representationDim();
  MatrixXd Jf = MatrixXd::Random(c, repDim);
  VectorXd x = Space.getZero().value();
  Space.retractation(x, x, VectorXd::Random(dim));
  MatrixXd expectedRes;
  expectedRes = Jf * Space.diffRetractation(x);
  MatrixXd J(c, dim);
  Space.applyDiffRetractation(J, Jf, x);
  BOOST_CHECK(expectedRes.isApprox(J));
}

BOOST_AUTO_TEST_CASE(RealApplyInvDiff)
{
  int c = 5;
  RealSpace Space(7);
  Index dim = Space.dim();
  Index repDim = Space.representationDim();
  MatrixXd Jf = MatrixXd::Random(c, dim);
  VectorXd x = Space.getZero().value();
  Space.retractation(x, x, VectorXd::Random(dim));
  MatrixXd expectedRes;
  expectedRes = Jf * Space.diffPseudoLog0(x);
  MatrixXd J(c, repDim);
  Space.applyDiffPseudoLog0(J, Jf, x);
  BOOST_CHECK(expectedRes.isApprox(J));
}

BOOST_AUTO_TEST_CASE(RealTransport)
{
  int c = 5;
  RealSpace Space(7);
  Index dim = Space.dim();
  MatrixXd H = MatrixXd::Random(dim, c);
  MatrixXd Hout(dim, c);
  VectorXd v = VectorXd::Random(dim);
  VectorXd x = Space.getZero().value();
  Space.retractation(x, x, v);
  MatrixXd expectedRes = H;
  Space.applyTransport(Hout, H, x, v);
  BOOST_CHECK(expectedRes.isApprox(Hout));
}

BOOST_AUTO_TEST_CASE(RealInvTransport)
{
  int c = 3;
  RealSpace Space(9);
  Index dim = Space.dim();
  MatrixXd H = MatrixXd::Random(dim, c);
  MatrixXd Hout(dim, c);
  VectorXd v = VectorXd::Random(dim);
  VectorXd x = Space.getZero().value();
  Space.retractation(x, x, v);
  MatrixXd expectedRes = H;
  Space.applyInvTransport(Hout, H, x, v);
  BOOST_CHECK(expectedRes.isApprox(Hout));
}

BOOST_AUTO_TEST_CASE(RealInvTransportOnTheRight)
{
  int c = 3;
  RealSpace Space(9);
  Index dim = Space.dim();
  MatrixXd H = MatrixXd::Random(c, dim);
  MatrixXd Hout(c, dim);
  VectorXd v = VectorXd::Random(dim);
  VectorXd x = Space.getZero().value();
  Space.retractation(x, x, v);
  MatrixXd expectedRes = H;
  Space.applyInvTransportOnTheRight(Hout, H, x, v);
  BOOST_CHECK(expectedRes.isApprox(Hout));
}

BOOST_AUTO_TEST_CASE(RealLimitMap)
{
  RealSpace Space(9);
  Index dim = Space.dim();
  VectorXd res(dim);
  Space.limitMap(res);
  VectorXd expectedRes(dim);
  double i = std::numeric_limits<double>::infinity();
  expectedRes << i, i, i, i, i, i, i, i, i;
  BOOST_CHECK_EQUAL(expectedRes, res);
}

BOOST_AUTO_TEST_CASE(isSameTopology)
{
  RealSpace R9(9);
  RealSpace R9b(9);
  RealSpace R3(3);
  S2 s2;
  SO3<ExpMapMatrix> so3;
  CartesianProduct cp(so3, R9);

  BOOST_CHECK(R9.isSameTopology(R9));
  BOOST_CHECK(R9.isSameTopology(R9b));
  BOOST_CHECK(!R9.isSameTopology(R3));
  BOOST_CHECK(!R9.isSameTopology(s2));
  BOOST_CHECK(!R9.isSameTopology(so3));
  BOOST_CHECK(!R9.isSameTopology(cp));
}

#if EIGEN_WORLD_VERSION > 3 ||                               \
    (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION > 2) || \
    (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION == 2 && \
     EIGEN_MINOR_VERSION > 0)
BOOST_AUTO_TEST_CASE(RealNoAllocation)
{
  // We only test here that the operations on the manifold do not create
  // temporary. Passing arguments that are not recognize by the Eigen::Ref will
  // create temporaries, but this is the user's fault.
  RealSpace R(4);
  Index c = 3;
  VectorXd x = VectorXd::Random(R.representationDim());
  VectorXd y = VectorXd::Random(R.representationDim());
  VectorXd v = VectorXd::Random(R.dim());
  VectorXd z(R.representationDim());
  MatrixXd J0 = MatrixXd::Random(c, R.representationDim());
  MatrixXd J1(c, R.representationDim());
  MatrixXd J2(c, R.representationDim());
  MatrixXd H0 = MatrixXd::Random(R.dim(), R.dim());
  MatrixXd H1 = MatrixXd::Random(R.dim(), R.dim());
  MatrixXd H2 = MatrixXd::Random(R.dim(), R.dim());

  internal::set_is_malloc_allowed(false);
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
  internal::set_is_malloc_allowed(true);
}
#endif

BOOST_AUTO_TEST_CASE(distance)
{
  RealSpace R9(9);
  VectorXd x(9), y(9);
  double res, expRes;

  x << 0.1, 0.2, 0.4, 1.5, 0.2, 2.44, 0, -1, -4.3;
  y << 0.1, 0.2, 0.4, 1.5, 0.2, 2.44, 0, -1, -4.3;
  expRes = sqrt((y-x).transpose()*(y-x));
  res = R9.distance(x,y);
  BOOST_CHECK_EQUAL(res, expRes);

  y << 1.2, 0.43, -9.3, 1.2, 1, 0, 0, -0.78, 0.01;
  expRes = sqrt((y-x).transpose()*(y-x));
  res = R9.distance(x,y);
  BOOST_CHECK_EQUAL(res, expRes);
}

BOOST_AUTO_TEST_CASE(squaredDistance)
{
  RealSpace R9(9);
  VectorXd x(9), y(9);
  double res, expRes;

  x << 0.1, 0.2, 0.4, 1.5, 0.2, 2.44, 0, -1, -4.3;
  y << 0.1, 0.2, 0.4, 1.5, 0.2, 2.44, 0, -1, -4.3;
  expRes = (y-x).transpose()*(y-x);
  res = R9.squaredDistance(x,y);
  BOOST_CHECK_EQUAL(res, expRes);

  y << 1.2, 0.43, -9.3, 1.2, 1, 0, 0, -0.78, 0.01;
  expRes = (y-x).transpose()*(y-x);
  res = R9.squaredDistance(x,y);
  BOOST_CHECK_EQUAL(res, expRes);
}
