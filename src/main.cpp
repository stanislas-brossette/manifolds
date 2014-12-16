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

#include <pgsolver/solver/ExampleGeometricProblem.h>
#include <pgsolver/solver/Problem.h>
#include <pgsolver/solver/Solver.h>

//#include <EigenQP/QLD.h>
//#include <EigenQP/QuadProg.h>
//#include <EigenQP/LSSOL.h>

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
    ExampleGeometricProblem myProb;
    Solver mySolver;
    Eigen::VectorXd v0(3);
    v0 << 1,2,3;
    Point x0 = myProb.M().createPoint(v0);
    mySolver.solve(myProb, x0);
  }
  {
  }
#ifdef _WIN32
  system("pause");
#endif
  return 0;
}

