#include <iostream>
#include <stdexcept>

#include <pgsolver/SO3.h>
#include <pgsolver/RealSpace.h>
#include <pgsolver/CartesianProduct.h>
#include <pgsolver/Point.h>
#include <pgsolver/ExpMapMatrix.h>

#define BOOST_TEST_MODULE PGSolver 

#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>

using boost::test_tools::output_test_stream;

using namespace pgs;

BOOST_AUTO_TEST_CASE(CartProdConstructor)
{
  RealSpace R3(3);
  RealSpace R2(2);
  SO3<ExpMapMatrix> RotSpace;
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
  SO3<ExpMapMatrix> RotSpace;
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
  SO3<ExpMapMatrix> RotSpace;
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
  SO3<ExpMapMatrix> RotSpace;
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
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2R3R2(R2, R3);
  R2R3R2.multiply(R2);
  CartesianProduct SO3R2R3R2(RotSpace, R2R3R2);
  Point x = SO3R2R3R2.getIdentity();
  Point y = SO3R2R3R2.getIdentity();
  Eigen::VectorXd vx(10);
  Eigen::VectorXd vy(10);
  vx << 1,0.1,1,1,2,3,4,5,6,7;
  vy << 0.07,3,0.01,4,6,2,1,4,6,4;
  x = x + vx;
  y = y + vy;
  Point z = x+(y-x);
  std::cout << "x = " << std::endl << x << std::endl;
  std::cout << "y = " << std::endl << y << std::endl;
  std::cout << "z = x+(y-x) = " << std::endl << z << std::endl;

  BOOST_CHECK_EQUAL(z.value().size(), 16);
  BOOST_CHECK_CLOSE(y.value()[0], z.value()[0], 1e-8);
  BOOST_CHECK_CLOSE(y.value()[1], z.value()[1], 1e-8);
  BOOST_CHECK_CLOSE(y.value()[2], z.value()[2], 1e-8);
  BOOST_CHECK_CLOSE(y.value()[3], z.value()[3], 1e-8);
  BOOST_CHECK_CLOSE(y.value()[4], z.value()[4], 1e-8);
  BOOST_CHECK_CLOSE(y.value()[5], z.value()[5], 1e-8);
  BOOST_CHECK_CLOSE(y.value()[6], z.value()[6], 1e-8);
  BOOST_CHECK_CLOSE(y.value()[7], z.value()[7], 1e-8);
  BOOST_CHECK_CLOSE(y.value()[8], z.value()[8], 1e-8);
  BOOST_CHECK_CLOSE(y.value()[9], z.value()[9], 1e-8); 
}

BOOST_AUTO_TEST_CASE(CartProInvMap)
{
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2SO3(R2, RotSpace);
  CartesianProduct R2R3(R2, R3);
  CartesianProduct R2SO3R2R3(R2SO3, R2R3);
  Point x = R2SO3R2R3.getIdentity();
  Point Id = R2SO3R2R3.getIdentity();
  Eigen::VectorXd vx(10);
  Eigen::VectorXd vy(10);
  vx << -7,2,1,0.1,1,3,4,5,6,7;
  vy << 4,6,0.07,3,0.01,2,1,4,6,4;
  x = x + vx +vy;
  Eigen::VectorXd x0 = x.invMap();
  std::cout << "x = " << std::endl << x << std::endl;
  std::cout << "x0 = " << std::endl << x0.transpose() << std::endl;

  Point newX = Id.increment(x0); 
  std::cout << "newX = Id.increment(x.InvMap()) =" << std::endl << newX << std::endl;

  BOOST_CHECK_EQUAL(newX.value().size(), 16);
  BOOST_CHECK_CLOSE(newX.value()[0], x.value()[0], 1e-8);
  BOOST_CHECK_CLOSE(newX.value()[1], x.value()[1], 1e-8);
  BOOST_CHECK_CLOSE(newX.value()[2], x.value()[2], 1e-8);
  BOOST_CHECK_CLOSE(newX.value()[3], x.value()[3], 1e-8);
  BOOST_CHECK_CLOSE(newX.value()[4], x.value()[4], 1e-8);
  BOOST_CHECK_CLOSE(newX.value()[5], x.value()[5], 1e-8);
  BOOST_CHECK_CLOSE(newX.value()[6], x.value()[6], 1e-8);
  BOOST_CHECK_CLOSE(newX.value()[7], x.value()[7], 1e-8);
  BOOST_CHECK_CLOSE(newX.value()[8], x.value()[8], 1e-8);
  BOOST_CHECK_CLOSE(newX.value()[9], x.value()[9], 1e-8); 
  BOOST_CHECK_CLOSE(newX.value()[10], x.value()[10], 1e-8);
  BOOST_CHECK_CLOSE(newX.value()[11], x.value()[11], 1e-8);
  BOOST_CHECK_CLOSE(newX.value()[12], x.value()[12], 1e-8);
  BOOST_CHECK_CLOSE(newX.value()[13], x.value()[13], 1e-8);
  BOOST_CHECK_CLOSE(newX.value()[14], x.value()[14], 1e-8);
  BOOST_CHECK_CLOSE(newX.value()[15], x.value()[15], 1e-8);
}
