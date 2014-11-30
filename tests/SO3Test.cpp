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
  vy << 1,2,3;
  x.increment(vy);
  x.increment(vy);
  BOOST_CHECK_EQUAL(x.value().size(), 9);
  BOOST_CHECK_EQUAL(x.value()[0], 0.40779157774427449);
  BOOST_CHECK_EQUAL(x.value()[1], 0.83844036685008172);
  BOOST_CHECK_EQUAL(x.value()[2], -0.36155743714814598);
  BOOST_CHECK_EQUAL(x.value()[3], -0.65622239077139688);
  BOOST_CHECK_EQUAL(x.value()[4], 0.54445505980328812);
  BOOST_CHECK_EQUAL(x.value()[5], 0.52243742372160695);
  BOOST_CHECK_EQUAL(x.value()[6], 0.63488440126617318);
  BOOST_CHECK_EQUAL(x.value()[7], 0.024216504514447346);
  BOOST_CHECK_EQUAL(x.value()[8], 0.77222752990164401);
}

BOOST_AUTO_TEST_CASE(SO3Addition)
{
  SO3 RotSpace;
  Point y = RotSpace.getIdentity();
  Eigen::Vector3d vy;
  vy << 1,2,3;
  y=(y+vy)+vy;
  BOOST_CHECK_EQUAL(y.value().size(), 9);
  BOOST_CHECK_EQUAL(y.value()[0], 0.40779157774427449);
  BOOST_CHECK_EQUAL(y.value()[1], 0.83844036685008172);
  BOOST_CHECK_EQUAL(y.value()[2], -0.36155743714814598);
  BOOST_CHECK_EQUAL(y.value()[3], -0.65622239077139688);
  BOOST_CHECK_EQUAL(y.value()[4], 0.54445505980328812);
  BOOST_CHECK_EQUAL(y.value()[5], 0.52243742372160695);
  BOOST_CHECK_EQUAL(y.value()[6], 0.63488440126617318);
  BOOST_CHECK_EQUAL(y.value()[7], 0.024216504514447346);
  BOOST_CHECK_EQUAL(y.value()[8], 0.77222752990164401);
}
