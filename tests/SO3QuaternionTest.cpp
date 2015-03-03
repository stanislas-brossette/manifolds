#include <iostream>
#include <stdexcept>

#include <manifolds/defs.h>
#include <manifolds/utils.h>
#include <manifolds/pgs_assert.h>
#include <manifolds/SO3.h>
#include <manifolds/Point.h>
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

//BOOST_AUTO_TEST_CASE(SO3Constructor)
//{
//  SO3<ExpMapQuaternion> S;
//  Point x = S.getIdentity();
//  Eigen::VectorXd v(4);
//  v << 1,0,0,0;
//  Point y = S.createPoint(v);
//  BOOST_CHECK_EQUAL(x.value().size(), 9);
//  BOOST_CHECK_EQUAL(x.value()[0], 1);
//  BOOST_CHECK_EQUAL(x.value()[1], 0);
//  BOOST_CHECK_EQUAL(x.value()[2], 0);
//  BOOST_CHECK_EQUAL(x.value()[3], 0);
//  BOOST_CHECK_EQUAL(x.value()[4], 1);
//  BOOST_CHECK_EQUAL(x.value()[5], 0);
//  BOOST_CHECK_EQUAL(x.value()[6], 0);
//  BOOST_CHECK_EQUAL(x.value()[7], 0);
//  BOOST_CHECK_EQUAL(x.value()[8], 1);
//  BOOST_CHECK_EQUAL(y.value().size(), 9);
//  BOOST_CHECK_EQUAL(y.value()[0], 1);
//  BOOST_CHECK_EQUAL(y.value()[1], 0);
//  BOOST_CHECK_EQUAL(y.value()[2], 0);
//  BOOST_CHECK_EQUAL(y.value()[3], 0);
//  BOOST_CHECK_EQUAL(y.value()[4], 1);
//  BOOST_CHECK_EQUAL(y.value()[5], 0);
//  BOOST_CHECK_EQUAL(y.value()[6], 0);
//  BOOST_CHECK_EQUAL(y.value()[7], 0);
//  BOOST_CHECK_EQUAL(y.value()[8], 1);
//}
//
//BOOST_AUTO_TEST_CASE(SO3Increment)
//{
//  SO3<ExpMapMatrix> S;
//  Eigen::VectorXd x = S.getIdentity().value();
//  Eigen::Vector3d vy;
//  vy << 0.1,0.2,0.3;
//  S.plus(x, x, vy);
//  S.plus(x, x, vy);
//  BOOST_CHECK_EQUAL(x.size(), 9);
//  BOOST_CHECK_CLOSE(x[0],  0.751909095300295, 1e-12);
//  BOOST_CHECK_CLOSE(x[1],  0.583715086608147, 1e-12);
//  BOOST_CHECK_CLOSE(x[2], -0.306446422838863, 1e-12);
//  BOOST_CHECK_CLOSE(x[3], -0.507379423623623, 1e-12);
//  BOOST_CHECK_CLOSE(x[4],  0.809160842538688, 1e-12);
//  BOOST_CHECK_CLOSE(x[5],  0.296352579515415, 1e-12);
//  BOOST_CHECK_CLOSE(x[6],  0.420949917315650, 1e-12);
//  BOOST_CHECK_CLOSE(x[7], -0.067345590561841, 1e-12);
//  BOOST_CHECK_CLOSE(x[8],  0.904580421269344, 1e-12);
//}
//
//BOOST_AUTO_TEST_CASE(SO3Addition)
//{
//  SO3<ExpMapMatrix> S;
//  Eigen::VectorXd y = S.getIdentity().value();
//  Eigen::Vector3d vy;
//  vy << 0.1,0.2,0.3;
//  S.plus(y, y, vy);
//  S.plus(y, y, vy);
//  BOOST_CHECK_EQUAL(y.size(), 9);
//  BOOST_CHECK_CLOSE(y[0],  0.751909095300295, 1e-12);
//  BOOST_CHECK_CLOSE(y[1],  0.583715086608147, 1e-12);
//  BOOST_CHECK_CLOSE(y[2], -0.306446422838863, 1e-12);
//  BOOST_CHECK_CLOSE(y[3], -0.507379423623623, 1e-12);
//  BOOST_CHECK_CLOSE(y[4],  0.809160842538688, 1e-12);
//  BOOST_CHECK_CLOSE(y[5],  0.296352579515415, 1e-12);
//  BOOST_CHECK_CLOSE(y[6],  0.420949917315650, 1e-12);
//  BOOST_CHECK_CLOSE(y[7], -0.067345590561841, 1e-12);
//  BOOST_CHECK_CLOSE(y[8],  0.904580421269344, 1e-12);
//}
//
//BOOST_AUTO_TEST_CASE(SO3Substraction)
//{
//  SO3<ExpMapMatrix> S;
//  Eigen::Vector3d v( 0.2, 0.4, 0.6);
//  Eigen::VectorXd R1 = S.getIdentity().value();
//  Eigen::VectorXd R2(9);
//  S.plus(R1, R1, v);
//  S.plus(R2, R1, v);
//  Eigen::Vector3d d;
//  S.minus(d,R2,R1);
//  BOOST_CHECK_CLOSE(d[0], 0.2, 1e-8);
//  BOOST_CHECK_CLOSE(d[1], 0.4, 1e-8);
//  BOOST_CHECK_CLOSE(d[2], 0.6, 1e-8);
//}
//
//BOOST_AUTO_TEST_CASE(SO3PointInvMap)
//{
//  SO3<ExpMapMatrix> S;
//  Eigen::VectorXd x = S.getIdentity().value();
//  Eigen::VectorXd vy = Eigen::VectorXd::Random(S.dim());
//  S.plus(x, x, vy);
//  Eigen::VectorXd z(S.dim());
//  S.invMap(z, x); 
//  BOOST_CHECK(z.isApprox(vy));
//}
//
//BOOST_AUTO_TEST_CASE(SO3Diff)
//{
//  SO3<ExpMapMatrix> S;
//  Eigen::MatrixXd J;
//  Eigen::MatrixXd Jtest(9,3);
//  Jtest <<         0,  0.003580526006716, -0.558259982176135,
//                   0,  0.646068748944272,  0.634553549147103,
//                   0, -0.763270824459509,  0.534497507539106,
//  -0.003580526006716,                  0, -0.829658346630838,
//  -0.646068748944272,                  0, -0.424189774632061,
//   0.763270824459509,                  0, -0.362946363755562,
//   0.558259982176135,  0.829658346630838,                  0,
//  -0.634553549147103,  0.424189774632061,                  0,
//  -0.534497507539106,  0.362946363755562,                  0;
//  Eigen::Vector3d v(0.680375, -0.211234, 0.566198);
//  Eigen::VectorXd x = S.getIdentity().value();
//  S.plus(x, x, v);
//  J = S.diffMap(x);
//  BOOST_CHECK(J.isApprox(Jtest));
//}
//
//BOOST_AUTO_TEST_CASE(SO3ApplyDiff)
//{
//  int c = 5;
//  SO3<ExpMapMatrix> S;
//  Index dim = S.dim();
//  Index repDim = S.representationDim();
//  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(c,repDim);
//  Eigen::VectorXd x = S.getIdentity().value();
//  S.plus(x, x, Eigen::VectorXd::Random(dim));
//  Eigen::MatrixXd expectedRes;
//  expectedRes = Jf*S.diffMap(x);
//  Eigen::MatrixXd J(c,dim);
//  S.applyDiffMap(J, Jf, x);
//  BOOST_CHECK(expectedRes.isApprox(J));
//}
//
//BOOST_AUTO_TEST_CASE(SO3invDiff)
//{
//  SO3<ExpMapMatrix> S;
//  Eigen::MatrixXd J;
//  Eigen::MatrixXd Jtest(3,9);
//  Jtest <<
//  -0.064043491813865,                 0,                  0,                  0, -0.064043491813865, 0.545030346992499,                 0, -0.545030346992499, -0.064043491813865,
//  -0.110117993664377,                 0, -0.545030346992499,                  0, -0.110117993664377,                 0, 0.545030346992499,                  0, -0.110117993664377,
//  -0.042109988599266, 0.545030346992499,                  0, -0.545030346992499, -0.042109988599266,                 0,                 0,                  0, -0.042109988599266;
//  Eigen::Vector3d v(0.3403857, 0.58526775, 0.223811);
//  Eigen::VectorXd x = S.getIdentity().value();
//  S.plus(x,x,v);
//  J = S.diffInvMap(x);
//  BOOST_CHECK(J.isApprox(Jtest));
//}
//
////BOOST_AUTO_TEST_CASE(SO3invDiffSmallValue)
////{
////  SO3<ExpMapMatrix> S;
////  Eigen::MatrixXd J;
////  Eigen::MatrixXd Jtest(3,9);
////  Jtest <<
////  -0.064043491813865,                 0,                  0,                  0, -0.064043491813865, 0.545030346992499,                 0, -0.545030346992499, -0.064043491813865,
////  -0.110117993664377,                 0, -0.545030346992499,                  0, -0.110117993664377,                 0, 0.545030346992499,                  0, -0.110117993664377,
////  -0.042109988599266, 0.545030346992499,                  0, -0.545030346992499, -0.042109988599266,                 0,                 0,                  0, -0.042109988599266;
////  Eigen::Vector3d v( 1.0e-08*0.081125768865785, 1.0e-08*0.929385970968730, 1.0e-08*0.775712678608402);
////  Point x = S.getIdentity();
////  x.increment(v);
////  J = S.diffInvMap(x.value());
////  std::cout << "J = " << J << std::endl; 
////  BOOST_CHECK(J.isApprox(Jtest));
////}
//
//BOOST_AUTO_TEST_CASE(SO3ApplyInvDiff)
//{
//  int c = 5;
//  SO3<ExpMapMatrix> S;
//  Index dim = S.dim();
//  Index repDim = S.representationDim();
//  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(c,dim);
//  Eigen::VectorXd x = S.getIdentity().value();
//  S.plus(x, x, Eigen::VectorXd::Random(dim));
//  Eigen::MatrixXd expectedRes;
//  expectedRes = Jf*S.diffInvMap(x);
//  Eigen::MatrixXd J(c,repDim);
//  S.applyDiffInvMap(J, Jf, x);
//  BOOST_CHECK(expectedRes.isApprox(J));
//}
//
//BOOST_AUTO_TEST_CASE(SO3Transport)
//{
//  int c = 4;
//  SO3<ExpMapMatrix> S;
//  Index dim = S.dim();
//  Eigen::MatrixXd H(dim,c);
//  H <<  1, 2, 3, 4,
//        5, 6, 7, 8,
//        9,10,11,12;
//  Eigen::MatrixXd Hout(dim,c);
//  Eigen::VectorXd v(dim);
//  v <<  0.083549465660115, 0.164064455761495, 0.287252050630289;
//  Eigen::VectorXd x = S.getIdentity().value();
//  S.plus(x, x, Eigen::VectorXd::Random(3));
//  Eigen::MatrixXd expectedRes(dim,c);
//  expectedRes << 1.126248257109656, 1.969921592423433, 2.813594927737210, 3.657268263050987,
//                 4.539510349826134, 5.725092676723538, 6.910675003620942, 8.096257330518345,
//                 9.226289104899047, 10.165762281434207, 11.105235457969370, 12.044708634504529;
//  S.applyTransport(Hout, H, x, v);
//  BOOST_CHECK(expectedRes.isApprox(Hout));
//}
//
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
//  Eigen::VectorXd x = S.getIdentity().value();
//  S.plus(x, x, Eigen::VectorXd::Random(3));
//  Eigen::MatrixXd expectedRes(r,dim);
//  expectedRes <<  0.667168954696934, 1.299987987788895,  3.444548855437121,
//                  2.972337006917136, 4.096292499301232,  7.168375023495865,
//                  5.277505059137337, 6.892597010813567, 10.892201191554610,
//                  7.582673111357540, 9.688901522325903, 14.616027359613355;
//  S.applyInvTransport(Hout, H, x, v);
//  BOOST_CHECK(expectedRes.isApprox(Hout));
//}
//
//#if   EIGEN_WORLD_VERSION > 3 \
//  || (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION > 2) \
//  || (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION == 2 && EIGEN_MINOR_VERSION > 0)
//BOOST_AUTO_TEST_CASE(SO3NoAllocation)
//{
//  //We only test here that the operations on the manifold do not create
//  //temporary. Passing arguments that are not recognize by the Eigen::Ref will
//  //create temporaries, but this is the user's fault.
//  const int r = 100;
//  SO3<ExpMapMatrix> S;
//  Index dim = S.dim();
//  Index repDim = S.representationDim();
//  Eigen::VectorXd x = Eigen::VectorXd::Random(repDim);
//  Eigen::VectorXd p = Eigen::VectorXd::Random(dim);
//  Eigen::VectorXd y = Eigen::VectorXd::Random(repDim);
//  Eigen::VectorXd z(repDim);
//  Eigen::VectorXd d(dim);
//  Eigen::MatrixXd J0 = Eigen::MatrixXd::Random(r, repDim);
//  Eigen::MatrixXd J1(r, dim);
//  Eigen::MatrixXd J2(r, repDim);
//  Eigen::MatrixXd H0 = Eigen::MatrixXd::Random(S.dim(), S.dim());
//  Eigen::MatrixXd H1 = Eigen::MatrixXd::Random(S.dim(), S.dim());
//  Eigen::MatrixXd H2 = Eigen::MatrixXd::Random(S.dim(), S.dim());
//
//  //The first call to the following methods might trigger a memory allocation
//  //depending on the size of the Ji and the initial buffer size inside S.
//  //However, subsequent calls should not require any allocation, what we check
//  //after.
//  S.applyDiffMap(J1, J0, x);
//  S.applyDiffInvMap(J2, J1, x);
//
//  Eigen::internal::set_is_malloc_allowed(false);
//  utils::set_is_malloc_allowed(false);
//  {
//    std::cout << "Memory allocation tests:" << std::endl;
//    S.plus(z, x, p);
//    std::cout << "- method 'plus' passed" << std::endl;
//    S.minus(d, x, y);
//    std::cout << "- method 'minus' passed" << std::endl;
//    S.invMap(d, x);
//    std::cout << "- method 'invMap' passed" << std::endl;
//    S.applyDiffMap(J1, J0, x);
//    std::cout << "- method 'applyDiffMap' passed" << std::endl;
//    S.applyDiffInvMap(J2, J1, x);
//    std::cout << "- method 'applyDiffInvMap' passed" << std::endl;
//    S.applyTransport(H1, H0, x, p);
//    std::cout << "- method 'applyTransport' passed" << std::endl;
//    S.applyInvTransport(H2, H0, x, p);
//    std::cout << "- method 'applyInvTransport' passed" << std::endl;
//  }
//  utils::set_is_malloc_allowed(true);
//  Eigen::internal::set_is_malloc_allowed(true);
//}
#endif

