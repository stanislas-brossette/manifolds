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
  vy << 1,2,3,4,5,6,7,0.1,0.2,0.3;
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
  BOOST_CHECK_CLOSE(x.value()[7],   0.751909095300295, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[8],   0.583715086608147, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[9],  -0.306446422838863, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[10], -0.507379423623623, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[11],  0.809160842538688, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[12],  0.296352579515415, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[13],  0.420949917315650, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[14], -0.067345590561841, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[15],  0.904580421269344, 1e-12);
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
  vy << 0.1,0.2,0.3,1,2,3,4,5,6,7;
  x = x + vy + vy;
  BOOST_CHECK_EQUAL(x.value().size(), 16);
  BOOST_CHECK_CLOSE(x.value()[0],  0.751909095300295, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[1],  0.583715086608147, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[2], -0.306446422838863, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[3], -0.507379423623623, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[4],  0.809160842538688, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[5],  0.296352579515415, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[6],  0.420949917315650, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[7], -0.067345590561841, 1e-12);
  BOOST_CHECK_CLOSE(x.value()[8],  0.904580421269344, 1e-12);
  BOOST_CHECK_EQUAL(x.value()[9], 2); 
  BOOST_CHECK_EQUAL(x.value()[10],4); 
  BOOST_CHECK_EQUAL(x.value()[11],6); 
  BOOST_CHECK_EQUAL(x.value()[12],8); 
  BOOST_CHECK_EQUAL(x.value()[13],10);
  BOOST_CHECK_EQUAL(x.value()[14],12);
  BOOST_CHECK_EQUAL(x.value()[15],14);
}

BOOST_AUTO_TEST_CASE(CartProSubstraction)
{
  RealSpace R3(3);
  RealSpace R2(2);
  SO3 RotSpace;
  CartesianProduct P(R2, R3);
  P.multiply(R2);
  CartesianProduct Q(RotSpace, P);
  Point xQ = Q.getIdentity();
  Eigen::VectorXd vxQ(10);
  vxQ << 0.1,0.2,0.3,1,2,3,4,5,6,7;
  Point yQ = xQ + vxQ;
  Eigen::VectorXd z = yQ-xQ;
  BOOST_CHECK_EQUAL(z.size(), 10);
  BOOST_CHECK_CLOSE(z[0], 0.1, 1e-12);
  BOOST_CHECK_CLOSE(z[1], 0.2, 1e-12);
  BOOST_CHECK_CLOSE(z[2], 0.3, 1e-12);
  BOOST_CHECK_CLOSE(z[3], 1,   1e-12);
  BOOST_CHECK_CLOSE(z[4], 2,   1e-12);
  BOOST_CHECK_CLOSE(z[5], 3,   1e-12);
  BOOST_CHECK_CLOSE(z[6], 4,   1e-12);
  BOOST_CHECK_CLOSE(z[7], 5,   1e-12);
  BOOST_CHECK_CLOSE(z[8], 6,   1e-12);
  BOOST_CHECK_CLOSE(z[9], 7,   1e-12); 
}
