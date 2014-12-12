#include <pgsolver/solver/ExampleGeometricProblem.h>

namespace pgs
{
  RealSpace ExampleGeometricProblem::R3 = RealSpace(3);   

  ExampleGeometricProblem::ExampleGeometricProblem()
    : Problem(R3)
  {
  }

  void ExampleGeometricProblem::getTangentLB(RefVec) const
  {
  }
  void ExampleGeometricProblem::getTangentUB(RefVec) const
  {
  }

  void ExampleGeometricProblem::evalObj(double& out) const
  {
    // x = [x1, x2, x3];
    // f(x) = x1^2+x2^2+x3^2;
    Eigen::Vector3d v = (x()+z()).value();
    out = v.squaredNorm();
  }
  void ExampleGeometricProblem::evalObjGrad(RefMat out) const
  {
    // x = [x1, x2, x3];
    // f(x) = x1^2+x2^2+x3^2;
    // df/dx(x) = [2.x1, 2.x2, 2.x3]^T;
    Eigen::Matrix<double,1,3> v = (x()+z()).value();
    out << 2*v[0], 2*v[1], 2*v[2]; 
    M().applyDiffMap(out, out, x().value());
  }

  void ExampleGeometricProblem::evalLinCstr(RefVec, size_t) const
  {
  }
  void ExampleGeometricProblem::evalLinCstrGrad(RefVec, size_t) const
  {
  }
  void ExampleGeometricProblem::getLinCstrLB(RefVec, size_t) const
  {
  }
  void ExampleGeometricProblem::getLinCstrUB(RefVec, size_t) const
  {
  }


  void ExampleGeometricProblem::evalNonLinCstr(RefVec, size_t) const
  {
  }
  void ExampleGeometricProblem::evalNonLinCstrGrad(RefVec, size_t) const
  {
  }
  void ExampleGeometricProblem::getNonLinCstrLB(RefVec, size_t) const
  {
  }
  void ExampleGeometricProblem::getNonLinCstrUB(RefVec, size_t) const
  {
  }

}
