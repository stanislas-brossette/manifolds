#include <iostream>
#include <stdexcept>
#include <pgsolver/RealSpace.h>
#include <pgsolver/SO3.h>
#include <pgsolver/CartesianProduct.h>
#include <pgsolver/Point.h>
#include <pgsolver/ExpMapMatrix.h>

using namespace pgs;

int main()
{
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2R3R2(R2, R3);
  R2R3R2.multiply(R2);
  CartesianProduct SO3R2R3R2(RotSpace, R2R3R2);
  //{
  //  Eigen::VectorXd v(3);
  //  v << 1,2,3;
  //  Point x = R3.createPoint(v);
  //  std::cout << x << std::endl;
  //}
  //{
  //  Eigen::VectorXd v(9);
  //  v << 1,0,0,0,1,0,0,0,1;
  //  Point x = RotSpace.createPoint(v);
  //  std::cout << x << std::endl;
  //}
  {
    Eigen::VectorXd v(16);
    v << 1,0,0,0,1,0,0,0,1,1,2,3,4,5,6,7;
    Point x = SO3R2R3R2.createPoint(v);
    std::cout << "x = " << std::endl << x << std::endl;
    Eigen::VectorXd vx(10);
    vx << 1,0,1,1,2,3,4,5,6,7;
    x.increment(vx);
    Point y = SO3R2R3R2.createPoint(x.value());

    //Point x = SO3R2R3R2.getIdentity();
    //Point y = SO3R2R3R2.getIdentity();
    //Eigen::VectorXd vx(10);
    //Eigen::VectorXd vy(10);
    //vx << 1,0,1,1,2,3,4,5,6,7;
    //vy << 0,3,0,4,6,2,1,4,6,4;
    //x = x + vx;
    //y = y + vy;
    //std::cout << "x = " << std::endl << x << std::endl;
    //std::cout << "y = " << std::endl << y << std::endl;
    //std::cout << "x+(y-x) = " << std::endl << x+(y-x) << std::endl;
  }
  return 0;
}

