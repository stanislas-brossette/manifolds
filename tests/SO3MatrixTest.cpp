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

#include <manifolds/defs.h>
#include <manifolds/utils.h>
#include <manifolds/mnf_assert.h>
#include <manifolds/RealSpace.h>
#include <manifolds/S2.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/SO3.h>
#include <manifolds/Point.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/ExpMapQuaternion.h>

#include <Eigen/Core>
#include <Eigen/LU>

#ifndef _WIN32
#define BOOST_TEST_MODULE Manifolds
#endif

#include <boost/test/unit_test.hpp>

using namespace mnf;
using namespace Eigen;
typedef Map<Matrix3d> toMat3;

BOOST_AUTO_TEST_CASE(RotSpaceConstructor)
{
  SO3<ExpMapMatrix> S;
  BOOST_CHECK_EQUAL(S.dim(), 3);
  BOOST_CHECK_EQUAL(S.representationDim(), 9);
  BOOST_CHECK_EQUAL(S.numberOfSubManifolds(), 1);
  BOOST_CHECK(S.isElementary());
}

BOOST_AUTO_TEST_CASE(SO3Constructor)
{
  SO3<ExpMapMatrix> S;
  Point x = S.getZero();
  VectorXd v(9);
  v << 1, 0, 0, 0, 1, 0, 0, 0, 1;
  Point y = S.createPoint(v);
  BOOST_CHECK_EQUAL(x.value().size(), 9);
  BOOST_CHECK_EQUAL(x.value()[0], 1);
  BOOST_CHECK_EQUAL(x.value()[1], 0);
  BOOST_CHECK_EQUAL(x.value()[2], 0);
  BOOST_CHECK_EQUAL(x.value()[3], 0);
  BOOST_CHECK_EQUAL(x.value()[4], 1);
  BOOST_CHECK_EQUAL(x.value()[5], 0);
  BOOST_CHECK_EQUAL(x.value()[6], 0);
  BOOST_CHECK_EQUAL(x.value()[7], 0);
  BOOST_CHECK_EQUAL(x.value()[8], 1);
  BOOST_CHECK_EQUAL(y.value().size(), 9);
  BOOST_CHECK_EQUAL(y.value()[0], 1);
  BOOST_CHECK_EQUAL(y.value()[1], 0);
  BOOST_CHECK_EQUAL(y.value()[2], 0);
  BOOST_CHECK_EQUAL(y.value()[3], 0);
  BOOST_CHECK_EQUAL(y.value()[4], 1);
  BOOST_CHECK_EQUAL(y.value()[5], 0);
  BOOST_CHECK_EQUAL(y.value()[6], 0);
  BOOST_CHECK_EQUAL(y.value()[7], 0);
  BOOST_CHECK_EQUAL(y.value()[8], 1);
}

BOOST_AUTO_TEST_CASE(RandomSO3PointConstructor)
{
  SO3<ExpMapMatrix> S;
  Point y = S.createRandomPoint();
  BOOST_CHECK(y.isInM());
}

BOOST_AUTO_TEST_CASE(RandomSO3IsInM)
{
  SO3<ExpMapMatrix> S;
  Point y = S.createRandomPoint();
  BOOST_CHECK(y.isInM());
  y.value() = VectorXd::Random(S.representationDim());
  BOOST_CHECK(!y.isInM());
}

BOOST_AUTO_TEST_CASE(testForceOnSo3)
{
  SO3<ExpMapMatrix> S;
  VectorXd rotValue(9);
  rotValue << 0.8788530623517696, 0.4728751994675595, -0.06329565959394243,
      -0.2643699005034564, 0.5931297977727257, 0.7604640679886711,
      0.3971471396129008, -0.6516027077311767, 0.6462879085784368;
  VectorXd randM(9);
  randM << 0.008323901360074014, 0.002714234559198019, 0.004345938588653662,
      -0.007167948892883933, 0.002139377525141173, -0.009673988567513409,
      -0.005142264587405261, -0.007255368464279626, 0.006083535084539808;
  Point rot = S.createPoint(rotValue);
  VectorXd perturbedRotVec;
  perturbedRotVec = rot.value() + randM;
  S.forceOnM(perturbedRotVec, perturbedRotVec);
  BOOST_CHECK(S.isInM(perturbedRotVec));
}

BOOST_AUTO_TEST_CASE(SO3LimitMap)
{
  SO3<ExpMapMatrix> S;
  Index dim = S.dim();
  VectorXd res(dim);
  S.limitMap(res);
  VectorXd expectedRes(dim);
  double i = 0.9 * M_PI / sqrt(3);
  expectedRes << i, i, i;
  BOOST_CHECK_EQUAL(expectedRes, res);
}

BOOST_AUTO_TEST_CASE(SO3Addition)
{
  SO3<ExpMapMatrix> S;
  VectorXd y = S.getZero().value();
  Vector3d vy;
  vy << 0.1, 0.2, 0.3;
  S.retractation(y, y, vy);
  S.retractation(y, y, vy);
  BOOST_CHECK_EQUAL(y.size(), 9);
  Matrix3d solution;
  solution << 0.751909095300295, -0.507379423623623, 0.420949917315650,
      0.583715086608147, 0.809160842538688, -0.067345590561841,
      -0.306446422838863, 0.296352579515415, 0.904580421269344;
  toMat3 yMat(y.data());
  BOOST_CHECK(yMat.isApprox(solution));
}

BOOST_AUTO_TEST_CASE(SO3Substraction)
{
  SO3<ExpMapMatrix> S;
  Vector3d v(0.2, 0.4, 0.6);
  VectorXd R1 = S.getZero().value();
  VectorXd R2(9);
  S.retractation(R1, R1, v);
  S.retractation(R2, R1, v);
  Vector3d d;
  S.pseudoLog(d, R1, R2);
  BOOST_CHECK_CLOSE(d[0], 0.2, 1e-8);
  BOOST_CHECK_CLOSE(d[1], 0.4, 1e-8);
  BOOST_CHECK_CLOSE(d[2], 0.6, 1e-8);
}

BOOST_AUTO_TEST_CASE(SO3PointpseudoLog0)
{
  SO3<ExpMapMatrix> S;
  VectorXd x = S.getZero().value();
  VectorXd vy = VectorXd::Random(S.dim());
  S.retractation(x, x, vy);
  VectorXd z(S.dim());
  S.pseudoLog0(z, x);
  BOOST_CHECK(z.isApprox(vy));
}

BOOST_AUTO_TEST_CASE(SO3Diff)
{
  SO3<ExpMapMatrix> S;
  MatrixXd J;
  MatrixXd Jtest(9, 3);
  Jtest << 0, 0.003580526006716, -0.558259982176135, 0, 0.646068748944272,
      0.634553549147103, 0, -0.763270824459509, 0.534497507539106,
      -0.003580526006716, 0, -0.829658346630838, -0.646068748944272, 0,
      -0.424189774632061, 0.763270824459509, 0, -0.362946363755562,
      0.558259982176135, 0.829658346630838, 0, -0.634553549147103,
      0.424189774632061, 0, -0.534497507539106, 0.362946363755562, 0;
  Vector3d v(0.680375, -0.211234, 0.566198);
  VectorXd x = S.getZero().value();
  S.retractation(x, x, v);
  J = S.diffRetractation(x);
  BOOST_CHECK(J.isApprox(Jtest));
}

BOOST_AUTO_TEST_CASE(SO3ApplyDiff)
{
  int c = 5;
  SO3<ExpMapMatrix> S;
  Index dim = S.dim();
  Index repDim = S.representationDim();
  MatrixXd Jf = MatrixXd::Random(c, repDim);
  VectorXd x = S.getZero().value();
  S.retractation(x, x, VectorXd::Random(dim));
  MatrixXd expectedRes;
  expectedRes = Jf * S.diffRetractation(x);
  MatrixXd J(c, dim);
  S.applyDiffRetractation(J, Jf, x);
  BOOST_CHECK(expectedRes.isApprox(J));
}

BOOST_AUTO_TEST_CASE(SO3invDiff)
{
  SO3<ExpMapMatrix> S;
  MatrixXd J;
  MatrixXd Jtest(3, 9);
  Jtest << -0.064043491813865, 0, 0, 0, -0.064043491813865, 0.545030346992499,
      0, -0.545030346992499, -0.064043491813865, -0.110117993664377, 0,
      -0.545030346992499, 0, -0.110117993664377, 0, 0.545030346992499, 0,
      -0.110117993664377, -0.042109988599266, 0.545030346992499, 0,
      -0.545030346992499, -0.042109988599266, 0, 0, 0, -0.042109988599266;
  Vector3d v(0.3403857, 0.58526775, 0.223811);
  VectorXd x = S.getZero().value();
  S.retractation(x, x, v);
  J = S.diffPseudoLog0(x);
  BOOST_CHECK(J.isApprox(Jtest));
}

// BOOST_AUTO_TEST_CASE(SO3invDiffSmallValue)
//{
//  SO3<ExpMapMatrix> S;
//  MatrixXd J;
//  MatrixXd Jtest(3,9);
//  Jtest <<
//  -0.064043491813865,                 0,                  0,
//  0, -0.064043491813865, 0.545030346992499,                 0,
//  -0.545030346992499, -0.064043491813865,
//  -0.110117993664377,                 0, -0.545030346992499,
//  0, -0.110117993664377,                 0, 0.545030346992499,
//  0, -0.110117993664377,
//  -0.042109988599266, 0.545030346992499,                  0,
//  -0.545030346992499, -0.042109988599266,                 0,
//  0,                  0, -0.042109988599266;
//  Vector3d v( 1.0e-08*0.081125768865785, 1.0e-08*0.929385970968730,
//  1.0e-08*0.775712678608402);
//  Point x = S.getZero();
//  x.increment(v);
//  J = S.diffPseudoLog0(x.value());
//  std::cout << "J = " << J << std::endl;
//  BOOST_CHECK(J.isApprox(Jtest));
//}

BOOST_AUTO_TEST_CASE(SO3ApplyInvDiff)
{
  int c = 5;
  SO3<ExpMapMatrix> S;
  Index dim = S.dim();
  Index repDim = S.representationDim();
  MatrixXd Jf = MatrixXd::Random(c, dim);
  VectorXd x = S.getZero().value();
  S.retractation(x, x, VectorXd::Random(dim));
  MatrixXd expectedRes;
  expectedRes = Jf * S.diffPseudoLog0(x);
  MatrixXd J(c, repDim);
  S.applyDiffPseudoLog0(J, Jf, x);
  BOOST_CHECK(expectedRes.isApprox(J));
}

BOOST_AUTO_TEST_CASE(isSameTopology)
{
  RealSpace R9(9);
  RealSpace R3(3);
  S2 s2;
  SO3<ExpMapMatrix> so3M;
  SO3<ExpMapQuaternion> so3Q;
  CartesianProduct cp(so3M, R9);

  BOOST_CHECK(!so3M.isSameTopology(R9));
  BOOST_CHECK(!so3M.isSameTopology(R3));
  BOOST_CHECK(!so3M.isSameTopology(s2));
  BOOST_CHECK(so3M.isSameTopology(so3M));
  BOOST_CHECK(!so3M.isSameTopology(so3Q));
  BOOST_CHECK(!so3M.isSameTopology(cp));
}

// BOOST_AUTO_TEST_CASE(SO3Transport)
//{
//  int c = 4;
//  SO3<ExpMapMatrix> S;
//  Index dim = S.dim();
//  MatrixXd H(dim,c);
//  H <<  1, 2, 3, 4,
//        5, 6, 7, 8,
//        9,10,11,12;
//  MatrixXd Hout(dim,c);
//  VectorXd v(dim);
//  v <<  0.083549465660115, 0.164064455761495, 0.287252050630289;
//  VectorXd x = S.getZero().value();
//  S.retractation(x, x, VectorXd::Random(3));
//  MatrixXd expectedRes(dim,c);
//  expectedRes << 1.126248257109656, 1.969921592423433, 2.813594927737210,
//  3.657268263050987,
//                 4.539510349826134, 5.725092676723538, 6.910675003620942,
//                 8.096257330518345,
//                 9.226289104899047, 10.165762281434207, 11.105235457969370,
//                 12.044708634504529;
//  S.applyTransport(Hout, H, x, v);
//  BOOST_CHECK(expectedRes.isApprox(Hout));
//}
//
// BOOST_AUTO_TEST_CASE(SO3InvTransport)
//{
//  int r = 4;
//  SO3<ExpMapMatrix> S;
//  Index dim = S.dim();
//  MatrixXd H(r,dim);
//  H <<  1, 2, 3,
//        4, 5, 6,
//        7, 8, 9,
//        10, 11, 12;
//  MatrixXd Hout(r,dim);
//  VectorXd v(dim);
//  v << 0.289466560559783, 0.047283924503264, 0.291177834528185;
//  VectorXd x = S.getZero().value();
//  S.retractation(x, x, VectorXd::Random(3));
//  MatrixXd expectedRes(r,dim);
//  expectedRes <<  0.667168954696934, 1.299987987788895,  3.444548855437121,
//                  2.972337006917136, 4.096292499301232,  7.168375023495865,
//                  5.277505059137337, 6.892597010813567, 10.892201191554610,
//                  7.582673111357540, 9.688901522325903, 14.616027359613355;
//  S.applyInvTransport(Hout, H, x, v);
//  BOOST_CHECK(expectedRes.isApprox(Hout));
//}

#if EIGEN_WORLD_VERSION > 3 ||                               \
    (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION > 2) || \
    (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION == 2 && \
     EIGEN_MINOR_VERSION > 0)
BOOST_AUTO_TEST_CASE(SO3NoAllocation)
{
  // We only test here that the operations on the manifold do not create
  // temporary. Passing arguments that are not recognize by the Ref will
  // create temporaries, but this is the user's fault.
  const int r = 100;
  SO3<ExpMapMatrix> S;
  Index dim = S.dim();
  Index repDim = S.representationDim();
  VectorXd x = VectorXd::Random(repDim);
  VectorXd p = VectorXd::Random(dim);
  VectorXd y = VectorXd::Random(repDim);
  VectorXd z(repDim);
  VectorXd d(dim);
  MatrixXd J0 = MatrixXd::Random(r, repDim);
  MatrixXd J1(r, dim);
  MatrixXd J2(r, repDim);
  MatrixXd H0 = MatrixXd::Random(S.dim(), S.dim());
  MatrixXd H1 = MatrixXd::Random(S.dim(), S.dim());
  MatrixXd H2 = MatrixXd::Random(S.dim(), S.dim());

  // The first call to the following methods might trigger a memory allocation
  // depending on the size of the Ji and the initial buffer size inside S.
  // However, subsequent calls should not require any allocation, what we check
  // after.
  S.applyDiffRetractation(J1, J0, x);
  S.applyDiffPseudoLog0(J2, J1, x);

  Eigen::internal::set_is_malloc_allowed(false);
  utils::set_is_malloc_allowed(false);
  {
    std::cout << "Memory allocation tests:" << std::endl;
    S.retractation(z, x, p);
    std::cout << "- method 'retractation' passed" << std::endl;
    S.pseudoLog(d, y, x);
    std::cout << "- method 'pseudoLog' passed" << std::endl;
    S.pseudoLog0(d, x);
    std::cout << "- method 'pseudoLog0' passed" << std::endl;
    S.applyDiffRetractation(J1, J0, x);
    std::cout << "- method 'applyDiffRetractation' passed" << std::endl;
    S.applyDiffPseudoLog0(J2, J1, x);
    std::cout << "- method 'applyDiffPseudoLog0' passed" << std::endl;
    S.applyTransport(H1, H0, x, p);
    std::cout << "- method 'applyTransport' passed" << std::endl;
    S.applyInvTransport(H2, H0, x, p);
    std::cout << "- method 'applyInvTransport' passed" << std::endl;
  }
  utils::set_is_malloc_allowed(true);
  Eigen::internal::set_is_malloc_allowed(true);
}
#endif

BOOST_AUTO_TEST_CASE(distance)
{
  SO3<ExpMapMatrix> so3;
  VectorXd x(9), y(9);
  toMat3 xMat(x.data());
  toMat3 yMat(y.data());
  double res, expRes;

  x = so3.getZero().value();
  y = so3.getZero().value();
  res = so3.distance(x, y);
  expRes = 0;
  BOOST_CHECK_CLOSE(res, expRes, 1e-9);

  Vector3d v(0.1, 0.5, -0.32);
  so3.retractation(x, x, v);
  res = so3.distance(x, y);
  expRes = v.norm();
  BOOST_CHECK_CLOSE(res, expRes, 1e-9);

  so3.retractation(y, y, v);
  v << -0.1, 0.42, 0;
  so3.retractation(x, x, v);
  res = so3.distance(x, y);
  expRes = v.norm();
  BOOST_CHECK_CLOSE(res, expRes, 1e-9);
}

BOOST_AUTO_TEST_CASE(squaredDistance)
{
  SO3<ExpMapMatrix> so3;
  VectorXd x(9), y(9);
  toMat3 xMat(x.data());
  toMat3 yMat(y.data());
  double res, expRes;

  x = so3.getZero().value();
  y = so3.getZero().value();
  res = so3.squaredDistance(x, y);
  expRes = 0;
  BOOST_CHECK_CLOSE(res, expRes, 1e-9);

  Vector3d v(0.1, 0.5, -0.32);
  so3.retractation(x, x, v);
  res = so3.squaredDistance(x, y);
  expRes = v.squaredNorm();
  BOOST_CHECK_CLOSE(res, expRes, 1e-9);

  so3.retractation(y, y, v);
  v << -0.1, 0.42, 0;
  so3.retractation(x, x, v);
  res = so3.squaredDistance(x, y);
  expRes = v.squaredNorm();
  BOOST_CHECK_CLOSE(res, expRes, 1e-9);
}
