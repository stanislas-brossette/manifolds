#include <iostream>
#include <stdexcept>

#include <pgsolver/SO3.h>
#include <pgsolver/RealSpace.h>
#include <pgsolver/CartesianProduct.h>
#include <pgsolver/Point.h>

#define BOOST_TEST_MODULE PGSolver 

#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>

using boost::test_tools::output_test_stream;

using namespace pgs;

BOOST_AUTO_TEST_CASE(CartProdConstructor)
{
  RealSpace R3(3);
  RealSpace R2(2);
  SO3 RotSpace;
  CartesianProduct P(R2, R3);
  P.multiply(R2);
  CartesianProduct Q(P, RotSpace);
  BOOST_CHECK_EQUAL(Q.dim(), 10);
  BOOST_CHECK_EQUAL(Q.representationDim(), 16);
  BOOST_CHECK_EQUAL(Q.numberOfSubmanifolds(), 2);
  BOOST_CHECK_EQUAL(P.numberOfSubmanifolds(), 3);
}

BOOST_AUTO_TEST_CASE(CartProdIdentity)
{
  RealSpace R3(3);
  RealSpace R2(2);
  SO3 RotSpace;
  CartesianProduct P(R2, R3);
  P.multiply(R2);
  CartesianProduct Q(P, RotSpace);
  Point x = Q.getIdentity();
  BOOST_CHECK_EQUAL(x.value().size(), 16);
  BOOST_CHECK_EQUAL(x.value()[0], 0);
  BOOST_CHECK_EQUAL(x.value()[1], 0);
  BOOST_CHECK_EQUAL(x.value()[2], 0);
  BOOST_CHECK_EQUAL(x.value()[3], 0);
  BOOST_CHECK_EQUAL(x.value()[4], 0);
  BOOST_CHECK_EQUAL(x.value()[5], 0);
  BOOST_CHECK_EQUAL(x.value()[6], 0);
  BOOST_CHECK_EQUAL(x.value()[7], 1);
  BOOST_CHECK_EQUAL(x.value()[8], 0);
  BOOST_CHECK_EQUAL(x.value()[9], 0);
  BOOST_CHECK_EQUAL(x.value()[10], 0);
  BOOST_CHECK_EQUAL(x.value()[11], 1);
  BOOST_CHECK_EQUAL(x.value()[12], 0);
  BOOST_CHECK_EQUAL(x.value()[13], 0);
  BOOST_CHECK_EQUAL(x.value()[14], 0);
  BOOST_CHECK_EQUAL(x.value()[15], 1);
}

BOOST_AUTO_TEST_CASE(CartProdIncrement)
{
  RealSpace R3(3);
  RealSpace R2(2);
  SO3 RotSpace;
  CartesianProduct P(R2, R3);
  P.multiply(R2);
  CartesianProduct Q(P, RotSpace);
  Point x = Q.getIdentity();
  Eigen::VectorXd vy(10);
  vy << 1,2,3,4,5,6,7,1,2,3;
  x.increment(vy);
  x.increment(vy);
  BOOST_CHECK_EQUAL(x.value().size(), 16);
  BOOST_CHECK_EQUAL(x.value()[0], 2);
  BOOST_CHECK_EQUAL(x.value()[1], 4);
  BOOST_CHECK_EQUAL(x.value()[2], 6);
  BOOST_CHECK_EQUAL(x.value()[3], 8);
  BOOST_CHECK_EQUAL(x.value()[4], 10);
  BOOST_CHECK_EQUAL(x.value()[5], 12);
  BOOST_CHECK_EQUAL(x.value()[6], 14);
  BOOST_CHECK_EQUAL(x.value()[7], 0.40779157774427449);
  BOOST_CHECK_EQUAL(x.value()[8], 0.83844036685008172);
  BOOST_CHECK_EQUAL(x.value()[9], -0.36155743714814598);
  BOOST_CHECK_EQUAL(x.value()[10], -0.65622239077139688);
  BOOST_CHECK_EQUAL(x.value()[11], 0.54445505980328812);
  BOOST_CHECK_EQUAL(x.value()[12], 0.52243742372160695);
  BOOST_CHECK_EQUAL(x.value()[13], 0.63488440126617318);
  BOOST_CHECK_EQUAL(x.value()[14], 0.024216504514447346);
  BOOST_CHECK_EQUAL(x.value()[15], 0.77222752990164401);
}

BOOST_AUTO_TEST_CASE(CartProdAddition)
{
  RealSpace R3(3);
  RealSpace R2(2);
  SO3 RotSpace;
  CartesianProduct P(R2, R3);
  P.multiply(R2);
  CartesianProduct Q(RotSpace, P);
  Point x = Q.getIdentity();
  Eigen::VectorXd vy(10);
  vy << 1,2,3,1,2,3,4,5,6,7;
  x = x + vy + vy;
  BOOST_CHECK_EQUAL(x.value().size(), 16);
  BOOST_CHECK_EQUAL(x.value()[0], 0.40779157774427449);
  BOOST_CHECK_EQUAL(x.value()[1], 0.83844036685008172);
  BOOST_CHECK_EQUAL(x.value()[2], -0.36155743714814598);
  BOOST_CHECK_EQUAL(x.value()[3], -0.65622239077139688);
  BOOST_CHECK_EQUAL(x.value()[4], 0.54445505980328812);
  BOOST_CHECK_EQUAL(x.value()[5], 0.52243742372160695);
  BOOST_CHECK_EQUAL(x.value()[6], 0.63488440126617318);
  BOOST_CHECK_EQUAL(x.value()[7], 0.024216504514447346);
  BOOST_CHECK_EQUAL(x.value()[8], 0.77222752990164401);
  BOOST_CHECK_EQUAL(x.value()[9], 2); 
  BOOST_CHECK_EQUAL(x.value()[10],4); 
  BOOST_CHECK_EQUAL(x.value()[11],6); 
  BOOST_CHECK_EQUAL(x.value()[12],8); 
  BOOST_CHECK_EQUAL(x.value()[13],10);
  BOOST_CHECK_EQUAL(x.value()[14],12);
  BOOST_CHECK_EQUAL(x.value()[15],14);
}
