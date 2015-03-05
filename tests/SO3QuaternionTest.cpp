#include <iostream>
#include <stdexcept>

#include <manifolds/defs.h>
#include <manifolds/utils.h>
#include <manifolds/pgs_assert.h>
#include <manifolds/SO3.h>
#include <manifolds/Point.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/ExpMapQuaternion.h>

#ifndef _WIN32
#define BOOST_TEST_MODULE PGSolver 
#endif

#include <boost/test/unit_test.hpp>

using namespace pgs;

BOOST_AUTO_TEST_CASE(RotSpaceConstructor)
{
  SO3<ExpMapQuaternion> S;
  BOOST_CHECK_EQUAL(S.dim(), 3);
  BOOST_CHECK_EQUAL(S.representationDim(), 4);
  BOOST_CHECK_EQUAL(S.numberOfSubmanifolds(), 1);
}

BOOST_AUTO_TEST_CASE(SO3Zero)
{
  SO3<ExpMapQuaternion> S;
  Point x = S.getZero();
  Eigen::Map<Eigen::Quaterniond> xQ(x.value().data());
  BOOST_CHECK(xQ.matrix().isApprox(Eigen::Matrix3d::Identity()));
}

BOOST_AUTO_TEST_CASE(SO3Constructor)
{
  SO3<ExpMapQuaternion> S;
  Point x = S.getZero();
  Eigen::VectorXd v(4);
  v << 0,0,0,1;
  Point y = S.createPoint(v);
  BOOST_CHECK_EQUAL(x.value().size(), 4);
  BOOST_CHECK_EQUAL(x.value()[0], 0);
  BOOST_CHECK_EQUAL(x.value()[1], 0);
  BOOST_CHECK_EQUAL(x.value()[2], 0);
  BOOST_CHECK_EQUAL(x.value()[3], 1);
  BOOST_CHECK_EQUAL(y.value().size(), 4);
  BOOST_CHECK_EQUAL(y.value()[0], 0);
  BOOST_CHECK_EQUAL(y.value()[1], 0);
  BOOST_CHECK_EQUAL(y.value()[2], 0);
  BOOST_CHECK_EQUAL(y.value()[3], 1);
}

BOOST_AUTO_TEST_CASE(SO3Addition)
{
  SO3<ExpMapQuaternion> S;
  Eigen::Vector4d x = S.getZero().value();
  Eigen::Map<Eigen::Quaterniond> xQ(x.data());

  Eigen::Vector3d vy;
  vy << 0.1,0.2,0.3;
  S.retractation(x, x, vy);
  S.retractation(x, x, vy);
  Eigen::Matrix3d solution;
  solution <<  0.751909095300295,-0.507379423623623, 0.420949917315650, 
               0.583715086608147, 0.809160842538688,-0.067345590561841,
              -0.306446422838863, 0.296352579515415, 0.904580421269344;
  BOOST_CHECK(xQ.matrix().isApprox(solution));
}

BOOST_AUTO_TEST_CASE(SO3InvMap)
{
  SO3<ExpMapQuaternion> S;
  Eigen::Vector4d x = S.getZero().value();

  Eigen::Vector3d v;
  v << 0.12364,-0.2234234,0.325843516;
  S.retractation(x, x, v);
  Eigen::Vector3d logX;
  S.invMap(logX, x);
  BOOST_CHECK(logX.isApprox(v));
}

BOOST_AUTO_TEST_CASE(SO3Substraction)
{
  SO3<ExpMapQuaternion> S;
  Eigen::Vector3d v( 0.268745647, -0.3544, 0.355546);
  Eigen::VectorXd q1 = S.getZero().value();
  Eigen::VectorXd q2(4);
  S.retractation(q1, q1, v);
  S.retractation(q2, q1, v);
  Eigen::Vector3d d;
  S.minus(d,q2,q1);
  BOOST_CHECK_CLOSE(d[0], v(0), 1e-8);
  BOOST_CHECK_CLOSE(d[1], v(1), 1e-8);
  BOOST_CHECK_CLOSE(d[2], v(2), 1e-8);
}

BOOST_AUTO_TEST_CASE(SO3Diff)
{
  double prec = 1e-9;
  SO3<ExpMapQuaternion> S;
  Eigen::Vector3d v(0.680375, -0.211234, 0.566198);
  Eigen::Vector4d q = S.getZero().value();
  S.retractation(q, q, v);
  Eigen::Matrix<double, 4, 3> J;
  Eigen::Vector4d dqdvx, dqdvy, dqdvz; 
  Eigen::Vector4d qpdx, qpdy, qpdz; 
  Eigen::Vector3d dvx, dvy, dvz;
  dvx << prec, 0, 0;
  dvy << 0, prec, 0;
  dvz << 0, 0, prec;
  S.retractation(qpdx,q,dvx);
  S.retractation(qpdy,q,dvy);
  S.retractation(qpdz,q,dvz);
  J.col(0) = (qpdx-q)/prec;
  J.col(1) = (qpdy-q)/prec;
  J.col(2) = (qpdz-q)/prec;

  Eigen::Matrix<double, 4, 3> diffM = S.diffRetractation(q);

  BOOST_CHECK(J.isApprox(diffM, 1e-6));
}

BOOST_AUTO_TEST_CASE(SO3ApplyDiff)
{
  int c = 5;
  SO3<ExpMapQuaternion> S;
  Index dim = S.dim();
  Index repDim = S.representationDim();
  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(c,repDim);
  Eigen::VectorXd x = S.getZero().value();
  S.retractation(x, x, Eigen::VectorXd::Random(dim));
  Eigen::MatrixXd expectedRes;
  expectedRes = Jf*S.diffRetractation(x);
  Eigen::MatrixXd J(c,dim);
  S.applyDiffRetractation(J, Jf, x);
  BOOST_CHECK(expectedRes.isApprox(J));
}

BOOST_AUTO_TEST_CASE(SO3InvDiff)
{
  double prec = 1e-6;
  SO3<ExpMapQuaternion> S;
  Eigen::Vector4d q = S.getZero().value();
  Eigen::Vector4d qpdx, qpdy, qpdz, qpdw;

  Eigen::Vector3d v;
  v << 0.12364,-0.2234234,0.325843516;
  S.retractation(q, q, v);
  qpdx = q;
  qpdx(0) += prec;
  qpdy = q;
  qpdy(1) += prec;
  qpdz = q;
  qpdz(2) += prec;
  qpdw = q;
  qpdw(3) += prec;
  Eigen::Vector3d logQ, logQpdx, logQpdy, logQpdz, logQpdw;
  S.invMap(logQ, q);
  S.invMap(logQpdx, qpdx);
  S.invMap(logQpdy, qpdy);
  S.invMap(logQpdz, qpdz);
  S.invMap(logQpdw, qpdw);

  Eigen::Matrix<double, 3, 4> J;
  J.col(0) = (logQpdx-logQ)/prec;
  J.col(1) = (logQpdy-logQ)/prec;
  J.col(2) = (logQpdz-logQ)/prec;
  J.col(3) = (logQpdw-logQ)/prec;

  Eigen::Matrix<double, 3, 4> invDiffM = S.diffInvMap(q);

  BOOST_CHECK(J.isApprox(invDiffM, 1e-6));
}
//BOOST_AUTO_TEST_CASE(SO3invDiffSmallValue)
//{
//  SO3<ExpMapMatrix> S;
//  Eigen::MatrixXd J;
//  Eigen::MatrixXd Jtest(3,9);
//  Jtest <<
//  -0.064043491813865,                 0,                  0,                  0, -0.064043491813865, 0.545030346992499,                 0, -0.545030346992499, -0.064043491813865,
//  -0.110117993664377,                 0, -0.545030346992499,                  0, -0.110117993664377,                 0, 0.545030346992499,                  0, -0.110117993664377,
//  -0.042109988599266, 0.545030346992499,                  0, -0.545030346992499, -0.042109988599266,                 0,                 0,                  0, -0.042109988599266;
//  Eigen::Vector3d v( 1.0e-08*0.081125768865785, 1.0e-08*0.929385970968730, 1.0e-08*0.775712678608402);
//  Point x = S.getZero();
//  x.increment(v);
//  J = S.diffInvMap(x.value());
//  std::cout << "J = " << J << std::endl; 
//  BOOST_CHECK(J.isApprox(Jtest));
//}

BOOST_AUTO_TEST_CASE(SO3ApplyInvDiff)
{
  int c = 5;
  SO3<ExpMapQuaternion> S;
  Index dim = S.dim();
  Index repDim = S.representationDim();
  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(c,dim);
  Eigen::VectorXd x = S.getZero().value();
  S.retractation(x, x, Eigen::VectorXd::Random(dim));
  Eigen::MatrixXd expectedRes;
  expectedRes = Jf*S.diffInvMap(x);
  Eigen::MatrixXd J(c,repDim);
  S.applyDiffInvMap(J, Jf, x);
  BOOST_CHECK(expectedRes.isApprox(J));
}

//BOOST_AUTO_TEST_CASE(SO3Transport)
//{
//  std::cout << "SO3Transport" << std::endl;
//  
//  int c = 4;
//  SO3<ExpMapQuaternion> S;
//  Index dim = S.dim();
//  Eigen::MatrixXd H(dim,c);
//  H <<  1, 2, 3, 4,
//        5, 6, 7, 8,
//        9,10,11,12;
//  Eigen::MatrixXd Hout(dim,c);
//  Eigen::VectorXd v(dim);
//  v <<  0.083549465660115, 0.164064455761495, 0.287252050630289;
//  Eigen::VectorXd x = S.getZero().value();
//  S.retractation(x, x, Eigen::VectorXd::Random(3));
//  Eigen::MatrixXd expectedRes(dim,c);
//  expectedRes << 1.126248257109656, 1.969921592423433, 2.813594927737210, 3.657268263050987,
//                 4.539510349826134, 5.725092676723538, 6.910675003620942, 8.096257330518345,
//                 9.226289104899047, 10.165762281434207, 11.105235457969370, 12.044708634504529;
//  std::cout << "COUCOU" << std::endl;
//  S.applyTransport(Hout, H, x, v);
//  BOOST_CHECK(expectedRes.isApprox(Hout));
//}

//BOOST_AUTO_TEST_CASE(SO3InvTransport)
//{
//  int r = 4;
//  SO3<ExpMapMatrix> S;
//  Index dim = S.dim();
//  Eigen::MatrixXd H(r,dim);
//  H <<  1, 2, 3,
//        4, 5, 6,
//        7, 8, 9,
//        10, 11, 12;
//  Eigen::MatrixXd Hout(r,dim);
//  Eigen::VectorXd v(dim);
//  v << 0.289466560559783, 0.047283924503264, 0.291177834528185;
//  Eigen::VectorXd x = S.getZero().value();
//  S.retractation(x, x, Eigen::VectorXd::Random(3));
//  Eigen::MatrixXd expectedRes(r,dim);
//  expectedRes <<  0.667168954696934, 1.299987987788895,  3.444548855437121,
//                  2.972337006917136, 4.096292499301232,  7.168375023495865,
//                  5.277505059137337, 6.892597010813567, 10.892201191554610,
//                  7.582673111357540, 9.688901522325903, 14.616027359613355;
//  S.applyInvTransport(Hout, H, x, v);
//  BOOST_CHECK(expectedRes.isApprox(Hout));
//}

Eigen::Vector3d rotatePointMatrix(Point& Rp, RefVec P)
{
  Eigen::Map<Eigen::Matrix3d> R(Rp.value().data());
  Eigen::Vector3d res;
  res = R*P;
  return res;
}
Eigen::Matrix<double,3,9> diffRotatePointMatrix(Point& , RefVec P)
{
  Eigen::Matrix<double,3,9> res;
  res << P(0), 0, 0, P(1), 0, 0, P(2), 0, 0, 
         0, P(0), 0, 0, P(1), 0, 0, P(2), 0,
         0, 0, P(0), 0, 0, P(1), 0, 0, P(2);
  return res;
}

Eigen::Vector3d rotatePointQuaternion(Point& Qp, RefVec P)
{
  Eigen::Map<Eigen::Quaterniond> Q(Qp.value().data());
  Eigen::Vector3d res;
  res <<  P.x() - 2*P.x()*Q.y()*Q.y() - 2*P.x()*Q.z()*Q.z() - 2*P.y()*Q.w()*Q.z() + 2*P.y()*Q.x()*Q.y() + 2*P.z()*Q.w()*Q.y() + 2*P.z()*Q.x()*Q.z(),
          P.y() - 2*P.y()*Q.x()*Q.x() - 2*P.y()*Q.z()*Q.z() + 2*P.x()*Q.w()*Q.z() + 2*P.x()*Q.x()*Q.y() - 2*P.z()*Q.w()*Q.x() + 2*P.z()*Q.y()*Q.z(),
          P.z() - 2*P.z()*Q.x()*Q.x() - 2*P.z()*Q.y()*Q.y() - 2*P.x()*Q.w()*Q.y() + 2*P.y()*Q.w()*Q.x() + 2*P.x()*Q.x()*Q.z() + 2*P.y()*Q.y()*Q.z();
  return res;
}

Eigen::Matrix<double,3,4> diffRotatePointQuaternion(Point& Qp, RefVec Pp)
{
  Eigen::Matrix<double,3,4> res;
  Eigen::Map<Eigen::Vector3d> P(Pp.data());
  Eigen::Map<Eigen::Quaterniond> Q(Qp.value().data());
  res.col(0) << 2*P.y()*Q.y() + 2*P.z()*Q.z(),
                2*P.x()*Q.y() - 4*P.y()*Q.x() - 2*P.z()*Q.w(),
                2*P.y()*Q.w() + 2*P.x()*Q.z() - 4*P.z()*Q.x();
  res.col(1) << 2*P.y()*Q.x() - 4*P.x()*Q.y() + 2*P.z()*Q.w(),
                2*P.x()*Q.x() + 2*P.z()*Q.z(),
                2*P.y()*Q.z() - 2*P.x()*Q.w() - 4*P.z()*Q.y();
  res.col(2) << 2*P.z()*Q.x() - 4*P.x()*Q.z() - 2*P.y()*Q.w(),
                2*P.x()*Q.w() - 4*P.y()*Q.z() + 2*P.z()*Q.y(),
                2*P.x()*Q.x() + 2*P.y()*Q.y();
  res.col(3) << 2*P.z()*Q.y() - 2*P.y()*Q.z(),
                2*P.x()*Q.z() - 2*P.z()*Q.x(),
                2*P.y()*Q.x() - 2*P.x()*Q.y();
  return res;
}

BOOST_AUTO_TEST_CASE(SO3CompareMatrixQuaternion)
{
  SO3<ExpMapMatrix> SO3_M;
  SO3<ExpMapQuaternion> SO3_Q;
  Point x_M = SO3_M.getZero();
  Point x_Q = SO3_Q.getZero();
  Eigen::Vector3d v = Eigen::Vector3d::Random();
  Eigen::Vector3d v2 = Eigen::Vector3d::Random();
  SO3_M.retractation(x_M.value(),x_M.value(),v);
  SO3_M.retractation(x_M.value(),x_M.value(),v2);
  SO3_Q.retractation(x_Q.value(),x_Q.value(),v);
  SO3_Q.retractation(x_Q.value(),x_Q.value(),v2);
  Eigen::Map<Eigen::Quaterniond> xQ(x_Q.value().data());
  Eigen::Map<Eigen::Matrix3d> xM(x_M.value().data());

  // Check that the expMaps for quaternion and 
  // rotation matrix are identical
  BOOST_CHECK(xQ.matrix().isApprox(xM));

  Eigen::Vector3d logX_Q, logX_M;
  SO3_Q.invMap(logX_Q,x_Q.value());
  SO3_M.invMap(logX_M,x_M.value());
  // Check that the logarithm for quaternion and 
  // rotation matrix are identical
  BOOST_CHECK(logX_Q.isApprox(logX_M));

  Eigen::Vector3d P0 = Eigen::Vector3d::Random();
  Eigen::Vector3d P0_Q = rotatePointQuaternion(x_Q, P0);
  Eigen::Vector3d P0_M = rotatePointMatrix(x_M, P0);

  // Check that the rotate function for quaternion and 
  // rotation matrix are identical
  BOOST_CHECK(P0_M.isApprox(P0_Q));

  Eigen::Matrix<double, 3, 3> J_M;
  Eigen::Matrix<double, 3, 9> J_M0 = diffRotatePointMatrix( x_M, P0);
  SO3_M.applyDiffRetractation(J_M, J_M0, x_M.value());

  Eigen::Matrix<double, 3, 3> J_Q;
  Eigen::Matrix<double, 3, 4> J_Q0 = diffRotatePointQuaternion( x_Q, P0);
  SO3_Q.applyDiffRetractation(J_Q, J_Q0, x_Q.value());

  // Check that the rotate function for quaternion and 
  // rotation matrix are identical
  BOOST_CHECK(J_M.isApprox(J_Q));
}


#if   EIGEN_WORLD_VERSION > 3 \
  || (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION > 2) \
  || (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION == 2 && EIGEN_MINOR_VERSION > 0)
BOOST_AUTO_TEST_CASE(SO3NoAllocation)
{
  //We only test here that the operations on the manifold do not create
  //temporary. Passing arguments that are not recognize by the Eigen::Ref will
  //create temporaries, but this is the user's fault.
  const int r = 100;
  SO3<ExpMapQuaternion> S;
  Index dim = S.dim();
  Index repDim = S.representationDim();
  Eigen::VectorXd x = Eigen::VectorXd::Random(repDim);
  Eigen::VectorXd p = Eigen::VectorXd::Random(dim);
  Eigen::VectorXd y = Eigen::VectorXd::Random(repDim);
  Eigen::VectorXd z(repDim);
  Eigen::VectorXd d(dim);
  Eigen::MatrixXd J0 = Eigen::MatrixXd::Random(r, repDim);
  Eigen::MatrixXd J1(r, dim);
  Eigen::MatrixXd J2(r, repDim);
  Eigen::MatrixXd H0 = Eigen::MatrixXd::Random(S.dim(), S.dim());
  Eigen::MatrixXd H1 = Eigen::MatrixXd::Random(S.dim(), S.dim());
  Eigen::MatrixXd H2 = Eigen::MatrixXd::Random(S.dim(), S.dim());

  //The first call to the following methods might trigger a memory allocation
  //depending on the size of the Ji and the initial buffer size inside S.
  //However, subsequent calls should not require any allocation, what we check
  //after.
  S.applyDiffRetractation(J1, J0, x);
  S.applyDiffInvMap(J2, J1, x);

  Eigen::internal::set_is_malloc_allowed(false);
  utils::set_is_malloc_allowed(false);
  {
    std::cout << "Memory allocation tests:" << std::endl;
    S.retractation(z, x, p);
    std::cout << "- method 'retractation' passed" << std::endl;
    S.minus(d, x, y);
    std::cout << "- method 'minus' passed" << std::endl;
    S.invMap(d, x);
    std::cout << "- method 'invMap' passed" << std::endl;
    S.applyDiffRetractation(J1, J0, x);
    std::cout << "- method 'applyDiffRetractation' passed" << std::endl;
    S.applyDiffInvMap(J2, J1, x);
    std::cout << "- method 'applyDiffInvMap' passed" << std::endl;
    S.applyTransport(H1, H0, x, p);
    std::cout << "- method 'applyTransport' passed" << std::endl;
    S.applyInvTransport(H2, H0, x, p);
    std::cout << "- method 'applyInvTransport' passed" << std::endl;
  }
  utils::set_is_malloc_allowed(true);
  Eigen::internal::set_is_malloc_allowed(true);
}
#endif

