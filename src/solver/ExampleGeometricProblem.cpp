#include <pgsolver/solver/ExampleGeometricProblem.h>

namespace pgs
{
      
  ExampleGeometricProblem::ExampleGeometricProblem()
    : Problem(*new RealSpace(3))      
  {
    std::cout << "Creating Example Geometric Problem" << std::endl;
    std::cout << "z =" << z() << std::endl;
    std::cout << "x =" << x() << std::endl;
  }

  void ExampleGeometricProblem::getTangentLB(RefVec) const
  {
  }
  void ExampleGeometricProblem::getTangentUB(RefVec) const
  {
  }

  void ExampleGeometricProblem::evalObj(RefVec) const
  {
  }
  void ExampleGeometricProblem::evalObjGrad(RefVec) const
  {
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
