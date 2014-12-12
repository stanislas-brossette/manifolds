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

  void ExampleProblem::getUB(RefVec) const
  {
  }
  void ExampleProblem::getLB(RefVec) const
  {
  }

  void ExampleProblem::getCstrLB(RefVec, size_t) const
  {
  }
  void ExampleProblem::getCstrUB(RefVec, size_t) const
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

  void ExampleProblem::evalNonLinCstr(RefVec, size_t) const
  {
  }
  void ExampleProblem::evalNonLinCstrGrad(RefVec, size_t) const
  {
  }
}
