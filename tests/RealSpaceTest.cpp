#include <iostream>
#include <stdexcept>

#include <pgsolver/pgs_assert.h>
#include <pgsolver/RealSpace.h>
#include <pgsolver/Point.h>

#define BOOST_TEST_MODULE PGSolver 

#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>

using boost::test_tools::output_test_stream;

using namespace pgs;

BOOST_AUTO_TEST_CASE(RealSpaceConstructor)
{
  RealSpace R5(5);
  BOOST_CHECK_EQUAL(R5.dim(), 5);
  BOOST_CHECK_EQUAL(R5.representationDim(), 5);
  BOOST_CHECK_EQUAL(R5.numberOfSubmanifolds(), 1);
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

BOOST_AUTO_TEST_CASE(RealSpaceIdentity)
{
  RealSpace R3(3);
  Point x = R3.getIdentity();
  for(long i = 0; i < x.value().size(); ++i)
  {
    BOOST_CHECK_EQUAL(x.value()[i], 0);
  }
}

BOOST_AUTO_TEST_CASE(RealPointIncrement)
{
  RealSpace R3(3);
  Point x = R3.getIdentity();
  Eigen::Vector3d vy;
  vy << 1,2,3;
  x.increment(vy);
  x.increment(vy);
  BOOST_CHECK_EQUAL(x.value().size(), 3);
  BOOST_CHECK_EQUAL(x.value()[0], 2);
  BOOST_CHECK_EQUAL(x.value()[1], 4);
  BOOST_CHECK_EQUAL(x.value()[2], 6);
}

BOOST_AUTO_TEST_CASE(RealPointAddition)
{
  RealSpace R3(3);
  Point y = R3.getIdentity();
  Eigen::Vector3d vy;
  vy << 1,2,3;
  y=y+vy+vy;
  BOOST_CHECK_EQUAL(y.value().size(), 3);
  BOOST_CHECK_EQUAL(y.value()[0], 2);
  BOOST_CHECK_EQUAL(y.value()[1], 4);
  BOOST_CHECK_EQUAL(y.value()[2], 6);
}

BOOST_AUTO_TEST_CASE(RealPointSubstraction)
{
  RealSpace R3(3);
  Point x = R3.getIdentity();
  Eigen::Vector3d vy;
  vy << 1,2,3;
  x = x + vy;
  Point y = x + vy + vy;
  Eigen::Vector3d z = y-x; 
  BOOST_CHECK_EQUAL(z[0], 2);
  BOOST_CHECK_EQUAL(z[1], 4);
  BOOST_CHECK_EQUAL(z[2], 6);
}

BOOST_AUTO_TEST_CASE(RealPointDiff)
{
  RealSpace R7(7);
  Eigen::MatrixXd J;
  Point x = R7.createPoint();
  J = R7.diffMap(x.value());
  bool test = J.isIdentity();
  BOOST_CHECK(test);
}

BOOST_AUTO_TEST_CASE(RealApplyDiff)
{
  RealSpace R7(7);
  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(5,7);
  Point x = R7.getIdentity();
  Eigen::MatrixXd expectedRes;
  expectedRes = Jf*R7.diffMap(x.value());
  Eigen::MatrixXd J(5,7);
  R7.applyDiffMap(J, Jf, x.value());
  bool test = expectedRes.isApprox(J);
  BOOST_CHECK(test);
}

BOOST_AUTO_TEST_CASE(RealApplyDiffGuaranteedResultTest)
{
  Index c = 3;
  RealSpace Space(7);
  Index dim = Space.dim();
  Index repDim = Space.representationDim();
  Eigen::MatrixXd Jf = Eigen::MatrixXd::Random(c,repDim);
  Eigen::MatrixXd Jres = Eigen::MatrixXd::Random(c,dim);
  Point x = Space.getIdentity();
  Space.applyDiffMap(Jres, Jf, x.value());
  
  bool worked = true;

  for (int i = 0; i<dim+1; ++i)
  {
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(c,repDim+dim);
    G.middleCols(dim,repDim) = Jf;
    Eigen::Map<Eigen::MatrixXd> Gf(G.data()+dim*c,c,repDim);
    Eigen::Map<Eigen::MatrixXd> Gres(G.data()+i*c,c,dim);
    Space.applyDiffMap(Gres,Gf,x.value());
    bool success = Jres.isApprox(Gres);
    worked = worked && success;
  }
  BOOST_CHECK(worked);
}
