#include <iostream>
#include <stdexcept>
#include <pgsolver/RealSpace.h>
#include <pgsolver/SO3.h>
#include <pgsolver/CartesianProduct.h>
#include <pgsolver/Point.h>

using namespace pgs;
int main()
{
  RealSpace R3(3);
  Eigen::Vector3d v;
  v << 1, 2, 3;
  Point x1 = R3.createPoint(v);
  Point x2 = R3.createPoint(4*v);
  std::cout << "x1=" <<  x1 << std::endl;
  std::cout << "x2=" <<  x2 << std::endl;
  std::cout << "x2-x1=" <<  x2-x1 << std::endl;

  //RealSpace R2(2);
  //SO3 RotSpace;
  //CartesianProduct P(R2, R3);
  //P.multiply(R2);
  //CartesianProduct Q(RotSpace, P);
  //Point x = Q.getIdentity();
  //Eigen::VectorXd vy(10);
  //vy << 1,2,3,1,2,3,4,5,6,7;
  //x = x + vy + vy;
  //std::cout << x;
  return 0;
}
