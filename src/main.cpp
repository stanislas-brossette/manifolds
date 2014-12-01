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
  std::cout << "x2-x1=" <<  (x2-x1).transpose() << std::endl;

  SO3 RotSpace;
  Point x = RotSpace.getIdentity();
  Eigen::Vector3d vy;
  vy << 0.1,0.2,0.3;
  x.increment(vy);
  x.increment(vy);

  v << 0.2, 0.4, 0.6;
  Point Rot1 = RotSpace.getIdentity();
  Rot1 = Rot1 + v;
  Point Rot2 = Rot1 + v;
  std::cout << "Rot1=" << std::endl <<  Rot1 << std::endl;
  std::cout << "Rot2=" << std::endl <<  Rot2 << std::endl;
  std::cout << "Rot2-Rot1=" << std::endl <<  Rot2-Rot1 << std::endl;
  
  RealSpace R2(2);
  CartesianProduct P(R2, R3);
  P.multiply(R2);
  CartesianProduct Q(RotSpace, P);
  Point xQ = Q.getIdentity();
  Eigen::VectorXd vxQ(10);
  vxQ << 0.1,0.2,0.3,1,2,3,4,5,6,7;
  Point yQ = xQ + vxQ + vxQ;
  std::cout << yQ-xQ;
  return 0;
}
