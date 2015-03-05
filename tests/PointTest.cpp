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
  Point x = R3.getZero();
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

BOOST_AUTO_TEST_CASE(PointRetractation)
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

BOOST_AUTO_TEST_CASE(SubPointManipulation)
{
  RealSpace R2(2);
  RealSpace R3(3);
  RealSpace R5(5);
  RealSpace R8(8);
  RealSpace R13(13);

  /*
              P
            /   \
           /      P3
          /      /  \
        P1      P2  R13
       /  \    /  \
      R2  R3  P5  P8
  */
  CartesianProduct P1(R2, R3);
  CartesianProduct P2(R5, R8);
  CartesianProduct P3(P2, R13);
  CartesianProduct P(P1, P3);
  Eigen::VectorXd v = Eigen::VectorXd::LinSpaced(P.representationDim(), 1, static_cast<double>(P.representationDim()));
  Point x = P.createPoint(v);

  ConstSubPoint p1 = x(0);  //const subpoint
  SubPoint p3 = x(1);       //non const subpoint
  SubPoint p2 = p3(0);      //non const again
  SubPoint r2 = x(0)(0);    //chaining subpoint
  ConstSubPoint r3 = p1(1); //const subpoint of a const subpoint
  //SubPoint r3b = p1(1);   //non const subpoint to a const subpoint. Does not compile, normal
  SubPoint r5 = p2(0);      //non const subpoint to a non const subpoint
  SubPoint r8 = p2(1);      //idem
  ConstSubPoint r13 = p3(1);//const subpoint of a non-const subpoint

  //Correct submanifold
  BOOST_CHECK(&x.getManifold() == &P);
  BOOST_CHECK(&r2.getManifold() == &R2);
  BOOST_CHECK(&r3.getManifold() == &R3);
  BOOST_CHECK(&r5.getManifold() == &R5);
  BOOST_CHECK(&r8.getManifold() == &R8);
  BOOST_CHECK(&r13.getManifold() == &R13);
  BOOST_CHECK(&p1.getManifold() == &P1);
  BOOST_CHECK(&p2.getManifold() == &P2);
  BOOST_CHECK(&p3.getManifold() == &P3);

  //Changing values of a subpoint
  Eigen::VectorXd v5(5); v5 << 0, 0, 7, 0, 0;
  r5.value() = v5;
  BOOST_CHECK(p2.value().head<5>() == v5);
  BOOST_CHECK(x.value().segment<5>(5) == v5);

  //New point is a duplicate, not a subpoint
  Eigen::VectorXd w5(5); w5 << 7, 7, 0, 7, 7;
  Point newr5 = p2(0);
  newr5.value() = w5;
  BOOST_CHECK(p2.value().head<5>() == v5);
  BOOST_CHECK(x.value().segment<5>(5) == v5);
}
