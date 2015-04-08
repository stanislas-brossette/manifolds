#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

#include <manifolds/defs.h>
#include <manifolds/utils.h>
#include <manifolds/pgs_assert.h>
#include <manifolds/S2.h>
#include <manifolds/Point.h>

#ifndef _WIN32
#define BOOST_TEST_MODULE PGSolver 
#endif

#include <boost/test/unit_test.hpp>

using namespace pgs;
using namespace Eigen;

BOOST_AUTO_TEST_CASE(S2Constructor)
{
  S2 S;
  BOOST_CHECK_EQUAL(S.dim(), 2);
  BOOST_CHECK_EQUAL(S.tangentDim(), 3);
  BOOST_CHECK_EQUAL(S.representationDim(), 3);
  BOOST_CHECK_EQUAL(S.numberOfSubmanifolds(), 1);
  BOOST_CHECK(S.isElementary());
}

BOOST_AUTO_TEST_CASE(S2PointConstructor)
{
  S2 S;
  Vector3d valX;
  S.rand(valX);
  Point x = S.createPoint(valX);
  Vector3d vy(1, 0, 0);
  Point y = S.createPoint(vy);
  BOOST_CHECK_EQUAL(x.value().size(), 3);
  BOOST_CHECK_EQUAL(y.value().size(), 3);
  BOOST_CHECK_EQUAL(y.value()[0], 1);
  BOOST_CHECK_EQUAL(y.value()[1], 0);
  BOOST_CHECK_EQUAL(y.value()[2], 0);
}

BOOST_AUTO_TEST_CASE(RandomS2Constructor)
{
  S2 Space;
  Point y = Space.createRandomPoint();
  BOOST_CHECK_EQUAL(y.value().size(), 3);
  BOOST_CHECK(y.isInM());
}

BOOST_AUTO_TEST_CASE(RandomS2IsInM)
{
  S2 S;
  Point y = S.createRandomPoint();
  BOOST_CHECK(y.isInM());
  y.value() = Eigen::VectorXd::Random(S.representationDim());
  BOOST_CHECK(!y.isInM());
}

BOOST_AUTO_TEST_CASE(testForceOnSo3)
{
  S2 S;
  Eigen::VectorXd perturbedRotVec(3);
  Eigen::VectorXd randM(3);
  Point rot = S.createRandomPoint();
  randM = 0.01*Eigen::VectorXd::Random(3);
  perturbedRotVec = rot.value() + randM;
  S.forceOnM(perturbedRotVec,perturbedRotVec);
  BOOST_CHECK(S.isInM(perturbedRotVec));
}

BOOST_AUTO_TEST_CASE(RandomTangentS2)
{
  S2 Space;
  Point x = Space.createRandomPoint();
  Vector3d v = Space.randVec(x.value());
  BOOST_CHECK(Space.isInTxM(x.value(), v));
}

BOOST_AUTO_TEST_CASE(S2Retractation)
{
  S2 Space;
  Vector3d vx;
  vx << 0, 0, 0;
  Vector3d valX0(1, 0, 0);
  Point x0 = Space.createPoint(valX0);
  Space.retractation(x0.value(),x0.value(),vx);
  BOOST_CHECK_EQUAL(x0.value().size(), 3);
  BOOST_CHECK(x0.value() == valX0);
  vx << 0, 1, 2;
  Space.retractation(x0.value(),x0.value(),vx);
  Vector3d xSol(0.408248290463863, 0.408248290463863, 0.816496580927726);
  BOOST_CHECK(x0.value().isApprox(xSol));
  vx << -1, -1, 1;
  Space.retractation(x0.value(),x0.value(),vx);
  xSol << -0.295875854768068, -0.295875854768068, 0.908248290463863;
  BOOST_CHECK(x0.value().isApprox(xSol));
}

BOOST_AUTO_TEST_CASE(S2Substraction)
{
  S2 Space;
  Vector3d x;
  Vector3d y;
  x << 1,0,0;
  y << 0,1,0;
  Vector3d z, zSol; 
  Space.pseudoLog(z,x,y);
  zSol << 0, M_PI/2, 0;
  BOOST_CHECK(z.isApprox(zSol));
}

BOOST_AUTO_TEST_CASE(S2ApplyDiff)
{
  int r = 5;
  S2 Space;
  Index tangentDim = Space.tangentDim();
  Index repDim = Space.representationDim();
  MatrixXd Jf = MatrixXd::Random(r,repDim);
  Jf.block(0, 0, 3, 3) << Matrix3d::Identity();
  Jf.block(3, 0, 2, 3) << 0, 0, 0,
                          1, 2, 3;
  Point x = Space.createPoint(Vector3d(1, 0, 0));
  MatrixXd expectedRes(r, tangentDim);
  expectedRes << 0, 0, 0,
                 0, 1, 0,
                 0, 0, 1,
                 0, 0, 0,
                 0, 2, 3;
  MatrixXd J(r, tangentDim);
  Space.applyDiffRetractation(J, Jf, x.value());
  BOOST_CHECK(expectedRes.isApprox(J));
  
  //then test random matrix
  Jf = MatrixXd::Random(r,repDim);
  x = Space.createRandomPoint();
  Space.applyDiffRetractation(J, Jf, x.value());
  
  for (Index i = 0; i < J.rows(); ++i)
    BOOST_CHECK(Space.isInTxM(J.row(i).transpose(),x.value()));
}

BOOST_AUTO_TEST_CASE(S2Transport)
{
  int c = 5;
  S2 Space;
  Index tangentDim = Space.tangentDim();
  MatrixXd H = MatrixXd::Random(tangentDim,c);
  MatrixXd Hout(tangentDim,c);
  Point x = Space.createRandomPoint();
  VectorXd v = Space.randVec(x.value());
  MatrixXd expectedRes = H;
  Space.applyTransport(Hout, H, x.value(), v);
  
  Eigen::Vector3d y;
  Space.retractation(y, x.value(), v);
  for (Index i = 0; i < Hout.cols(); ++i)
    BOOST_CHECK(Space.isInTxM(Hout.col(i).transpose(),y));

  H = MatrixXd::Zero(tangentDim,c);
  Space.applyTransport(Hout, H, x.value(), v);
  BOOST_CHECK_EQUAL(Hout, H);
}

BOOST_AUTO_TEST_CASE(S2InvTransport)
{
  int r = 7;
  S2 Space;
  Index tangentDim = Space.tangentDim();
  MatrixXd H = MatrixXd::Random(r, tangentDim);
  MatrixXd Hout(r, tangentDim);
  Point x = Space.createRandomPoint();
  VectorXd v = Space.randVec(x.value());
  MatrixXd expectedRes = H;
  Space.applyInvTransport(Hout, H, x.value(), v);
  
  for (Index i = 0; i < Hout.rows(); ++i)
    BOOST_CHECK(Space.isInTxM(Hout.row(i).transpose(),x.value()));

  H = MatrixXd::Zero(r, tangentDim);
  Space.applyInvTransport(Hout, H, x.value(), v);
  BOOST_CHECK(Hout == H);
}

BOOST_AUTO_TEST_CASE(S2LimitMap)
{
  S2 Space;
  Index tangentDim = Space.tangentDim();
  VectorXd res(tangentDim);
  Space.limitMap(res);
  VectorXd expectedRes(tangentDim);
  double i = std::numeric_limits<double>::infinity();
  expectedRes << i,i,i;
  BOOST_CHECK_EQUAL(expectedRes, res);
}

#if   EIGEN_WORLD_VERSION > 3 \
  || (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION > 2) \
  || (EIGEN_WORLD_VERSION == 3 && EIGEN_MAJOR_VERSION == 2 && EIGEN_MINOR_VERSION > 0)
BOOST_AUTO_TEST_CASE(S2NoAllocation)
{
  //We only test here that the operations on the manifold do not create
  //temporary. Passing arguments that are not recognize by the Eigen::Ref will
  //create temporaries, but this is the user's fault.
  S2 Space;
  Index c = 3;
  //VectorXd x = VectorXd::Random(Space.representationDim());
  //VectorXd y = VectorXd::Random(Space.representationDim());
  VectorXd x = Space.createRandomPoint().value();
  VectorXd y = Space.createRandomPoint().value();
  VectorXd vx = Space.randVec(x);
  VectorXd z(Space.representationDim());
  MatrixXd J0 = MatrixXd::Random(c, Space.representationDim());
  MatrixXd J1(c, Space.representationDim());
  MatrixXd J2(c, Space.representationDim());
  MatrixXd H0 = MatrixXd::Random(Space.tangentDim(), Space.tangentDim());
  MatrixXd H1 = MatrixXd::Random(Space.tangentDim(), Space.tangentDim());
  MatrixXd H2 = MatrixXd::Random(Space.tangentDim(), Space.tangentDim());

  internal::set_is_malloc_allowed(false);
  utils::set_is_malloc_allowed(false);
  {
    Space.retractation(z, x, vx);
    Space.pseudoLog(z, x, y);
    Space.applyDiffRetractation(J1, J0, x);
    Space.applyTransport(H1, H0, x, vx);
    Space.applyInvTransport(H2, H0, x, vx);
  }
  utils::set_is_malloc_allowed(true);
  internal::set_is_malloc_allowed(true);
}
#endif

