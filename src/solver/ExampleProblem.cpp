#include <pgsolver/solver/ExampleProblem.h>

namespace pgs
{
      
  ExampleProblem::ExampleProblem(Manifold& M)
    : Problem(M)
  {
  }
  ExampleProblem::ExampleProblem(Manifold& M, const Point& x)
    : Problem(M, x)
  {
  }

  void ExampleProblem::getTangentLB(RefVec) const
  {
  }
  void ExampleProblem::getTangentUB(RefVec) const
  {
  }

  void ExampleProblem::evalObj(RefVec) const
  {
  }
  void ExampleProblem::evalObjGrad(RefVec) const
  {
  }

  void ExampleProblem::evalLinCstr(RefVec, size_t) const
  {
  }
  void ExampleProblem::evalLinCstrGrad(RefVec, size_t) const
  {
  }
  void ExampleProblem::getLinCstrLB(RefVec, size_t) const
  {
  }
  void ExampleProblem::getLinCstrUB(RefVec, size_t) const
  {
  }


  void ExampleProblem::evalNonLinCstr(RefVec, size_t) const
  {
  }
  void ExampleProblem::evalNonLinCstrGrad(RefVec, size_t) const
  {
  }
  void ExampleProblem::getNonLinCstrLB(RefVec, size_t) const
  {
  }
  void ExampleProblem::getNonLinCstrUB(RefVec, size_t) const
  {
  }
}
