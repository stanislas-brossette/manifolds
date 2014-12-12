#include <iostream>
#include <stdexcept>
#define _USE_MATH_DEFINES
#include <math.h>

#include <Eigen/Core>

#include <pgsolver/manifolds/pgs_assert.h>
#include <pgsolver/manifolds/RealSpace.h>
#include <pgsolver/manifolds/SO3.h>
#include <pgsolver/manifolds/CartesianProduct.h>
#include <pgsolver/manifolds/Point.h>
#include <pgsolver/manifolds/ExpMapMatrix.h>
#include <pgsolver/manifolds/ReusableTemporaryMap.h>

#include <pgsolver/solver/Problem.h>
#include <pgsolver/solver/ExampleProblem.h>

using namespace pgs;

int main()
{
  std::cout << "Using: Eigen" << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION <<"." << EIGEN_MINOR_VERSION<< std::endl;
  RealSpace R2(2);
  RealSpace R3(3);
  SO3<ExpMapMatrix> RotSpace;
  CartesianProduct R2R3R2(R2, R3);
  R2R3R2.multiply(R2);
  CartesianProduct SO3R2R3R2(RotSpace, R2R3R2);
  CartesianProduct SO3R2R3R2SO3(SO3R2R3R2, RotSpace);
  CartesianProduct R3SO3(R3, RotSpace);
  CartesianProduct R3SO3R3SO3(R3, RotSpace);
  R3SO3R3SO3.multiply(R3);
  R3SO3R3SO3.multiply(RotSpace);
  {
    ExampleProblem myProb(R3SO3);
    Eigen::VectorXd v0(12);
    v0 << 3,4,5,1,0,0,0,1,0,0,0,1;
    Eigen::VectorXd v1(6);
    v1 << 3,4,5,1,0,0;
    Point x0 = R3SO3.createPoint(v0);
    myProb.setX(x0);
    std::cout << myProb.x()<<std::endl;
    myProb.setZ(v1);
    x0 = x0+v1;
    myProb.setX(x0);
    std::cout << myProb.x()<<std::endl;
  }
#ifdef _WIN32
  system("pause");
#endif
  return 0;
}

