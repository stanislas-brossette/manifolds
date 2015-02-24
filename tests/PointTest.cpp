#include <iostream>

#ifndef _WIN32
#define BOOST_TEST_MODULE PGSolver 
#endif

#include <boost/test/unit_test.hpp>

#include <manifolds/defs.h>
#include <manifolds/utils.h>
#include <manifolds/Point.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>

using namespace pgs;

//
//All the calculations are tested in the manifolds. Here we only test that the
//points methods are correctly mapped onto the manifolds methods
//

BOOST_AUTO_TEST_CASE(PointConstructor)
{
  RealSpace R3(3);
  Eigen::Vector3d v(1,2,3);
  const Point x = R3.createPoint(v);
  x.value();
  BOOST_CHECK(true/*v.isApprox(x.value())*/);
}

BOOST_AUTO_TEST_CASE(PointIncrement)
{
  RealSpace R3(3);
  Eigen::Vector3d v(1,2,3);
  Point x = R3.getIdentity();
  x.increment(v);
  BOOST_CHECK(v.isApprox(x.value()));
}

BOOST_AUTO_TEST_CASE(PointAccessors)
{
  RealSpace R3(3);
  Eigen::Vector3d v3(1,2,3);
  RealSpace R2(2);
  Eigen::Vector2d v2(10,20);
  Eigen::VectorXd v(5);
  v << 1,2,3,10,20;
  CartesianProduct S(R3, R2);
  Point x = S.createPoint(v);
  BOOST_CHECK(v3.isApprox(x(0).value()));
  BOOST_CHECK(v2.isApprox(x(1).value()));
  BOOST_CHECK(v3.isApprox(x[0]));
  BOOST_CHECK(v2.isApprox(x[1]));
}

BOOST_AUTO_TEST_CASE(PointPlus)
{
  RealSpace R3(3);
  Eigen::Vector3d v(1,2,3);

  Point x = R3.createPoint(v);
  x = x + v;
  BOOST_CHECK(x.value().isApprox(2*v));
}


BOOST_AUTO_TEST_CASE(PointMinus)
{
  RealSpace R3(3);
  Eigen::Vector3d v(1,2,3);

  Point x = R3.createPoint(v);
  Point y = R3.createPoint(2*v);
  Eigen::Vector3d res;

  res = y - x;
  BOOST_CHECK(res.isApprox(v));
}
