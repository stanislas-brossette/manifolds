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
  RealSpace R2(2);
  SO3 RotSpace;
  CartesianProduct P(R2, R3);
  P.multiply(R2);
  CartesianProduct Q(P, RotSpace);
  
  Eigen::VectorXd val = Eigen::VectorXd::Ones(7);
  Eigen::VectorXd p = Eigen::VectorXd::LinSpaced(7, 0, 6);
  Point x = P.createPoint();
  Point y = P.getIdentity();
  Point R = Q.getIdentity();
  try
  {
    P.multiply(R3);
  }
  catch (std::runtime_error e)
  {
    std::cout << e.what() << std::endl;
  }
  Point z = P.createPoint(val);
  x[0] = Eigen::Vector2d(7, 6);
  x[1] = Eigen::Vector3d(5, 4, 3);
  x[2] = Eigen::Vector2d(2, 1);
  std::cout << "p = " << p.transpose() <<std::endl;
  std::cout << "x = " << std::endl << x << std::endl;
  std::cout << "y = " << std::endl << y << std::endl;
  y = y+p;
  std::cout << "y = y+p = " << std::endl<< y << std::endl;
  std::cout << "R = " << std::endl<< R << std::endl;
  return 0;
}
