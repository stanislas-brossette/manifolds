#include <iostream>
#include <stdexcept>
#include <pgsolver/SO3.h>
#include <pgsolver/Point.h>

#define BOOST_TEST_MODULE PGSolver 

#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>

using boost::test_tools::output_test_stream;

using namespace pgs;

BOOST_AUTO_TEST_CASE(RotSpaceConstructor)
{
  SO3 RotSpace;
  BOOST_CHECK_EQUAL(RotSpace.dim(), 3);
  BOOST_CHECK_EQUAL(RotSpace.representationDim(), 9);
  BOOST_CHECK_EQUAL(RotSpace.numberOfSubmanifolds(), 1);
}

BOOST_AUTO_TEST_CASE(SO3Constructor)
{
  SO3 RotSpace;
  Point x = RotSpace.getIdentity();
  Eigen::VectorXd v(9);
  v << 1,2,3,4,5,6,7,8,9;
  Point y = RotSpace.createPoint(v);
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
  BOOST_CHECK_EQUAL(y.value()[1], 2);
  BOOST_CHECK_EQUAL(y.value()[2], 3);
  BOOST_CHECK_EQUAL(y.value()[3], 4);
  BOOST_CHECK_EQUAL(y.value()[4], 5);
  BOOST_CHECK_EQUAL(y.value()[5], 6);
  BOOST_CHECK_EQUAL(y.value()[6], 7);
  BOOST_CHECK_EQUAL(y.value()[7], 8);
  BOOST_CHECK_EQUAL(y.value()[8], 9);
}

BOOST_AUTO_TEST_CASE(SO3Increment)
{
  SO3 RotSpace;
  Point x = RotSpace.getIdentity();
  Eigen::Vector3d vy;
  vy << 0.1,0.2,0.3;
  x.increment(vy);
  x.increment(vy);
  BOOST_CHECK_EQUAL(x.value().size(), 9);
  BOOST_CHECK_EQUAL(x.value()[0], 0.75190909530029459);
  BOOST_CHECK_EQUAL(x.value()[1], 0.5837150866081473);
  BOOST_CHECK_EQUAL(x.value()[2], -0.30644642283886314);
  BOOST_CHECK_EQUAL(x.value()[3], -0.50737942362362254);
  BOOST_CHECK_EQUAL(x.value()[4], 0.80916084253868814);
  BOOST_CHECK_EQUAL(x.value()[5], 0.29635257951541538);
  BOOST_CHECK_EQUAL(x.value()[6], 0.42094991731565012);
  BOOST_CHECK_EQUAL(x.value()[7], -0.067345590561841279);
  BOOST_CHECK_EQUAL(x.value()[8], 0.90458042126934401);
}

BOOST_AUTO_TEST_CASE(SO3Addition)
{
  SO3 RotSpace;
  Point y = RotSpace.getIdentity();
  Eigen::Vector3d vy;
  vy << 0.1,0.2,0.3;
  y=(y+vy)+vy;
  BOOST_CHECK_EQUAL(y.value().size(), 9);
  BOOST_CHECK_EQUAL(y.value()[0], 0.75190909530029459);
  BOOST_CHECK_EQUAL(y.value()[1], 0.5837150866081473);
  BOOST_CHECK_EQUAL(y.value()[2], -0.30644642283886314);
  BOOST_CHECK_EQUAL(y.value()[3], -0.50737942362362254);
  BOOST_CHECK_EQUAL(y.value()[4], 0.80916084253868814);
  BOOST_CHECK_EQUAL(y.value()[5], 0.29635257951541538);
  BOOST_CHECK_EQUAL(y.value()[6], 0.42094991731565012);
  BOOST_CHECK_EQUAL(y.value()[7], -0.067345590561841279);
  BOOST_CHECK_EQUAL(y.value()[8], 0.90458042126934401);
}

BOOST_AUTO_TEST_CASE(SO3Substraction)
{
  SO3 RotSpace;
  Eigen::Vector3d v( 0.2, 0.4, 0.6);
  Point R1 = RotSpace.getIdentity();
  R1 = R1 + v;
  Point R2 = R1 + v;
  Eigen::Vector3d d = R2-R1;
  BOOST_CHECK_CLOSE(d[0], 0.2, 1e-8);
  BOOST_CHECK_CLOSE(d[1], 0.4, 1e-8);
  BOOST_CHECK_CLOSE(d[2], 0.6, 1e-8);
}
