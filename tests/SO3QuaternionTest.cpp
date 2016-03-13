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
#include <manifolds/Point.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/ExpMapQuaternion.h>

#ifndef _WIN32
#define BOOST_TEST_MODULE Manifolds
#endif

#include <boost/test/unit_test.hpp>

using namespace mnf;
using namespace Eigen;
typedef utils::ReverseQuaternion toQuat;
typedef utils::ConstReverseQuaternion toConstQuat;

BOOST_AUTO_TEST_CASE(RotSpaceConstructor)
{
  SO3<ExpMapQuaternion> S;
  BOOST_CHECK_EQUAL(S.dim(), 3);
  BOOST_CHECK_EQUAL(S.representationDim(), 4);
  BOOST_CHECK_EQUAL(S.numberOfSubManifolds(), 1);
  BOOST_CHECK(S.isElementary());
}

BOOST_AUTO_TEST_CASE(SO3Zero)
{
  SO3<ExpMapQuaternion> S;
  Point x = S.getZero();
  toQuat xQ(x.value().data());

  BOOST_CHECK(xQ.matrix().isApprox(Matrix3d::Identity()));
}

BOOST_AUTO_TEST_CASE(SO3Constructor)
{
  SO3<ExpMapQuaternion> S;
  Point x = S.getZero();
  VectorXd v(4);
  v << 1, 0, 0, 0;
  Point y = S.createPoint(v);
  BOOST_CHECK_EQUAL(x.value().size(), 4);
  BOOST_CHECK_EQUAL(x.value()[0], 1);
  BOOST_CHECK_EQUAL(x.value()[1], 0);
  BOOST_CHECK_EQUAL(x.value()[2], 0);
  BOOST_CHECK_EQUAL(x.value()[3], 0);
  BOOST_CHECK_EQUAL(y.value().size(), 4);
  BOOST_CHECK_EQUAL(y.value()[0], 1);
  BOOST_CHECK_EQUAL(y.value()[1], 0);
  BOOST_CHECK_EQUAL(y.value()[2], 0);
  BOOST_CHECK_EQUAL(y.value()[3], 0);
}

BOOST_AUTO_TEST_CASE(RandomSO3PointConstructor)
{
  SO3<ExpMapQuaternion> S;
  Point y = S.createRandomPoint();
  BOOST_CHECK_EQUAL(y.value().size(), 4);
  BOOST_CHECK(y.isInM());
}

BOOST_AUTO_TEST_CASE(RandomSO3IsInM)
{
  SO3<ExpMapQuaternion> S;
  Point y = S.createRandomPoint();
  BOOST_CHECK(y.isInM());
  y.value() = VectorXd::Random(S.representationDim());
  BOOST_CHECK(!y.isInM());
}

BOOST_AUTO_TEST_CASE(testForceOnSo3)
{
  // std::cout.precision(16);
  SO3<ExpMapQuaternion> S;
  VectorXd rotValue(4);
  VectorXd perturbedRotVec(4);
  VectorXd randM(4);
  rotValue << -0.02246981965305121, 0.1281118914339608, -0.134419124980123,
      0.9823512352094593;
  randM << 0.0002680182039123102, 0.009044594503494256, 0.008323901360074014,
      0.002714234559198019;
  // Point rot = S.createRandomPoint();
  Point rot = S.createPoint(rotValue);
  // std::cout << "rot.value() = \n" << rot.value() << std::endl;
  // randM = 0.01*VectorXd::Random(4);
  // std::cout << "randM = \n" << randM << std::endl;
  perturbedRotVec = rot.value() + randM;
  S.forceOnM(perturbedRotVec, perturbedRotVec);

  BOOST_CHECK(S.isInM(perturbedRotVec));
}

BOOST_AUTO_TEST_CASE(SO3LimitMap)
{
  SO3<ExpMapQuaternion> S;
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
  SO3<ExpMapQuaternion> S;
  Vector4d x = S.getZero().value();

  Vector3d vy;
  vy << 0.1, 0.2, 0.3;
  S.retractation(x, x, vy);
  S.retractation(x, x, vy);

  Matrix3d solution;
  solution << 0.751909095300295, -0.507379423623623, 0.420949917315650,
      0.583715086608147, 0.809160842538688, -0.067345590561841,
      -0.306446422838863, 0.296352579515415, 0.904580421269344;
  toQuat xQ(x.data());
  BOOST_CHECK(xQ.matrix().isApprox(solution));
}

BOOST_AUTO_TEST_CASE(SO3pseudoLog0)
{
  SO3<ExpMapQuaternion> S;
  Vector4d x = S.getZero().value();

  Vector3d v;
  v << 0.12364, -0.2234234, 0.325843516;
  S.retractation(x, x, v);
  Vector3d logX;
  S.pseudoLog0(logX, x);
  BOOST_CHECK(logX.isApprox(v));
}

BOOST_AUTO_TEST_CASE(SO3Substraction)
{
  SO3<ExpMapQuaternion> S;
  Vector3d v(0.268745647, -0.3544, 0.355546);
  VectorXd q1 = S.getZero().value();
  VectorXd q2(4);
  S.retractation(q1, q1, v);
  S.retractation(q2, q1, v);
  Vector3d d;
  S.pseudoLog(d, q1, q2);
  BOOST_CHECK_CLOSE(d[0], v(0), 1e-8);
  BOOST_CHECK_CLOSE(d[1], v(1), 1e-8);
  BOOST_CHECK_CLOSE(d[2], v(2), 1e-8);
}

BOOST_AUTO_TEST_CASE(SO3Diff)
{
  double prec = 1e-9;
  SO3<ExpMapQuaternion> S;
  Vector3d v(0.680375, -0.211234, 0.566198);
  Vector4d q = S.getZero().value();
  S.retractation(q, q, v);
  Matrix<double, 4, 3> J;
  Vector4d dqdvx, dqdvy, dqdvz;
  Vector4d qpdx, qpdy, qpdz;
  Vector3d dvx, dvy, dvz;
  dvx << prec, 0, 0;
  dvy << 0, prec, 0;
  dvz << 0, 0, prec;
  S.retractation(qpdx, q, dvx);
  S.retractation(qpdy, q, dvy);
  S.retractation(qpdz, q, dvz);
  J.col(0) = (qpdx - q) / prec;
  J.col(1) = (qpdy - q) / prec;
  J.col(2) = (qpdz - q) / prec;

  Matrix<double, 4, 3> diffM = S.diffRetractation(q);

  BOOST_CHECK(J.isApprox(diffM, 1e-6));
}

BOOST_AUTO_TEST_CASE(SO3ApplyDiff)
{
  int c = 5;
  SO3<ExpMapQuaternion> S;
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

BOOST_AUTO_TEST_CASE(SO3InvDiff)
{
  double prec = 1e-6;
  SO3<ExpMapQuaternion> S;
  Vector4d q = S.getZero().value();
  Vector4d qpdx, qpdy, qpdz, qpdw;

  Vector3d v;
  v << 0.12364, -0.2234234, 0.325843516;
  S.retractation(q, q, v);
  qpdx = q;
  qpdx(0) += prec;
  qpdy = q;
  qpdy(1) += prec;
  qpdz = q;
  qpdz(2) += prec;
  qpdw = q;
  qpdw(3) += prec;
  Vector3d logQ, logQpdx, logQpdy, logQpdz, logQpdw;
  S.pseudoLog0(logQ, q);
  S.pseudoLog0(logQpdx, qpdx);
  S.pseudoLog0(logQpdy, qpdy);
  S.pseudoLog0(logQpdz, qpdz);
  S.pseudoLog0(logQpdw, qpdw);

  Matrix<double, 3, 4> J;
  J.col(0) = (logQpdx - logQ) / prec;
  J.col(1) = (logQpdy - logQ) / prec;
  J.col(2) = (logQpdz - logQ) / prec;
  J.col(3) = (logQpdw - logQ) / prec;

  Matrix<double, 3, 4> invDiffM = S.diffPseudoLog0(q);

  BOOST_CHECK(J.isApprox(invDiffM, 1e-6));
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
  SO3<ExpMapQuaternion> S;
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

// BOOST_AUTO_TEST_CASE(SO3Transport)
//{
//  std::cout << "SO3Transport" << std::endl;
//
//  int c = 4;
//  SO3<ExpMapQuaternion> S;
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
//  std::cout << "COUCOU" << std::endl;
//  S.applyTransport(Hout, H, x, v);
//  BOOST_CHECK(expectedRes.isApprox(Hout));
//}

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

Vector3d rotatePointMatrix(Point& Rp, RefVec P)
{
  Map<Matrix3d> R(Rp.value().data());
  Vector3d res;
  res = R * P;
  return res;
}
Matrix<double, 3, 9> diffRotatePointMatrix(Point&, RefVec P)
{
  Matrix<double, 3, 9> res;
  res << P(0), 0, 0, P(1), 0, 0, P(2), 0, 0, 0, P(0), 0, 0, P(1), 0, 0, P(2), 0,
      0, 0, P(0), 0, 0, P(1), 0, 0, P(2);
  return res;
}

Vector3d rotatePointQuaternion(Point& Qp, RefVec P)
{
  toConstQuat Q(Qp.value().data());
  Vector3d res;
  res << P.x() - 2 * P.x() * Q.y() * Q.y() - 2 * P.x() * Q.z() * Q.z() -
             2 * P.y() * Q.w() * Q.z() + 2 * P.y() * Q.x() * Q.y() +
             2 * P.z() * Q.w() * Q.y() + 2 * P.z() * Q.x() * Q.z(),
      P.y() - 2 * P.y() * Q.x() * Q.x() - 2 * P.y() * Q.z() * Q.z() +
          2 * P.x() * Q.w() * Q.z() + 2 * P.x() * Q.x() * Q.y() -
          2 * P.z() * Q.w() * Q.x() + 2 * P.z() * Q.y() * Q.z(),
      P.z() - 2 * P.z() * Q.x() * Q.x() - 2 * P.z() * Q.y() * Q.y() -
          2 * P.x() * Q.w() * Q.y() + 2 * P.y() * Q.w() * Q.x() +
          2 * P.x() * Q.x() * Q.z() + 2 * P.y() * Q.y() * Q.z();
  return res;
}

Matrix<double, 3, 4> diffRotatePointQuaternion(Point& Qp, RefVec Pp)
{
  Matrix<double, 3, 4> res;
  Map<Vector3d> P(Pp.data());
  toConstQuat Q(Qp.value().data());

  res.col(0) << 2 * P.z() * Q.y() - 2 * P.y() * Q.z(),
      2 * P.x() * Q.z() - 2 * P.z() * Q.x(),
      2 * P.y() * Q.x() - 2 * P.x() * Q.y();
  res.col(1) << 2 * P.y() * Q.y() + 2 * P.z() * Q.z(),
      2 * P.x() * Q.y() - 4 * P.y() * Q.x() - 2 * P.z() * Q.w(),
      2 * P.y() * Q.w() + 2 * P.x() * Q.z() - 4 * P.z() * Q.x();
  res.col(2) << 2 * P.y() * Q.x() - 4 * P.x() * Q.y() + 2 * P.z() * Q.w(),
      2 * P.x() * Q.x() + 2 * P.z() * Q.z(),
      2 * P.y() * Q.z() - 2 * P.x() * Q.w() - 4 * P.z() * Q.y();
  res.col(3) << 2 * P.z() * Q.x() - 4 * P.x() * Q.z() - 2 * P.y() * Q.w(),
      2 * P.x() * Q.w() - 4 * P.y() * Q.z() + 2 * P.z() * Q.y(),
      2 * P.x() * Q.x() + 2 * P.y() * Q.y();
  return res;
}

BOOST_AUTO_TEST_CASE(SO3CompareMatrixQuaternion)
{
  SO3<ExpMapMatrix> SO3_M;
  SO3<ExpMapQuaternion> SO3_Q;
  Point x_M = SO3_M.getZero();
  Point x_Q = SO3_Q.getZero();
  Vector3d v = Vector3d::Random();
  Vector3d v2 = Vector3d::Random();
  SO3_M.retractation(x_M.value(), x_M.value(), v);
  SO3_M.retractation(x_M.value(), x_M.value(), v2);
  SO3_Q.retractation(x_Q.value(), x_Q.value(), v);
  SO3_Q.retractation(x_Q.value(), x_Q.value(), v2);
  toQuat xQ(x_Q.value().data());
  Map<Matrix3d> xM(x_M.value().data());

  // Check that the expMaps for quaternion and
  // rotation matrix are identical
  BOOST_CHECK(xQ.matrix().isApprox(xM));

  Vector3d logX_Q, logX_M;
  SO3_Q.pseudoLog0(logX_Q, x_Q.value());
  SO3_M.pseudoLog0(logX_M, x_M.value());
  // Check that the logarithm for quaternion and
  // rotation matrix are identical
  BOOST_CHECK(logX_Q.isApprox(logX_M));

  Vector3d P0 = Vector3d::Random();
  Vector3d P0_Q = rotatePointQuaternion(x_Q, P0);
  Vector3d P0_M = rotatePointMatrix(x_M, P0);

  // Check that the rotate function for quaternion and
  // rotation matrix are identical
  BOOST_CHECK(P0_M.isApprox(P0_Q));

  Matrix<double, 3, 3> J_M;
  Matrix<double, 3, 9> J_M0 = diffRotatePointMatrix(x_M, P0);
  SO3_M.applyDiffRetractation(J_M, J_M0, x_M.value());

  Matrix<double, 3, 3> J_Q;
  Matrix<double, 3, 4> J_Q0 = diffRotatePointQuaternion(x_Q, P0);
  SO3_Q.applyDiffRetractation(J_Q, J_Q0, x_Q.value());

  // Check that the rotate function for quaternion and
  // rotation matrix are identical
  BOOST_CHECK(J_M.isApprox(J_Q));
}

BOOST_AUTO_TEST_CASE(isSameTopology)
{
  RealSpace R9(9);
  RealSpace R3(3);
  S2 s2;
  SO3<ExpMapMatrix> so3M;
  SO3<ExpMapQuaternion> so3Q;
  CartesianProduct cp(so3M, R9);

  BOOST_CHECK(!so3Q.isSameTopology(R9));
  BOOST_CHECK(!so3Q.isSameTopology(R3));
  BOOST_CHECK(!so3Q.isSameTopology(s2));
  BOOST_CHECK(!so3Q.isSameTopology(so3M));
  BOOST_CHECK(so3Q.isSameTopology(so3Q));
  BOOST_CHECK(!so3Q.isSameTopology(cp));
}

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
  SO3<ExpMapQuaternion> S;
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
  SO3<ExpMapQuaternion> so3;
  VectorXd x(4), y(4);
  toQuat xQ(x.data());
  toQuat yQ(y.data());
  double res, expRes;

  so3.createRandomPoint(x);
  y = x;

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
  SO3<ExpMapQuaternion> so3;
  VectorXd x(4), y(4);
  toQuat xQ(x.data());
  toQuat yQ(y.data());
  double res, expRes;

  so3.createRandomPoint(x);
  y = x;

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

BOOST_AUTO_TEST_CASE(diffDistance)
{
  SO3<ExpMapQuaternion> M;
  VectorXd x(4), y(4), xTy(4);
  M.createRandomPoint(x);
  M.createRandomPoint(y);
  toConstQuat xQ(x.data());
  toConstQuat yQ(y.data());
  toQuat xTyQ(xTy.data());
  Eigen::Matrix<double, 1, 4> res, expRes;
  Eigen::Matrix<double, 3, 4> resLog, expResLog;

  xTyQ = xQ.inverse() * yQ;  // TODO double-check that formula
  xTyQ.writeChanges();

  double deltaFD = 1e-10;
  double deltaRes = 1e-3;

  resLog = M.diffPseudoLog0(xTy);
  expResLog = utils::FDLogarithm(M, xTy, deltaFD);
  BOOST_CHECK(resLog.isApprox(expResLog, deltaRes));

  res = M.derivDistanceX(x, y);
  expRes = utils::FDDerivDistanceX(M, x, y, deltaFD);
  BOOST_CHECK(res.isApprox(expRes, deltaRes));

  res = M.derivDistanceY(x, y);
  expRes = utils::FDDerivDistanceY(M, x, y, deltaFD);
  BOOST_CHECK(res.isApprox(expRes, deltaRes));

  res = M.derivSquaredDistanceX(x, y);
  expRes = utils::FDDerivSquaredDistanceX(M, x, y, deltaFD);
  BOOST_CHECK(res.isApprox(expRes, deltaRes));

  res = M.derivSquaredDistanceY(x, y);
  expRes = utils::FDDerivSquaredDistanceY(M, x, y, deltaFD);
  BOOST_CHECK(res.isApprox(expRes, deltaRes));
}
