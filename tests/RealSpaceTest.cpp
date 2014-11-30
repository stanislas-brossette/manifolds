#include <iostream>
#include <stdexcept>
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
    BOOST_CHECK_EQUAL(x.value()[0], 0);
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