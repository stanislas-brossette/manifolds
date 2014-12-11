#include <pgsolver/solver/ExampleProblem.h>

namespace pgs
{
      
  ExampleProblem::ExampleProblem(const Manifold& M)
    : Problem(M)
  {
  }


  void ExampleProblem::getUB(RefVec)
  {
  }
  void ExampleProblem::getLB(RefVec)
  {
  }

  void ExampleProblem::getCstrLB(RefVec)
  {
  }
  void ExampleProblem::getCstrUB(RefVec)
  {
  }

  void ExampleProblem::evalObj(RefVec)
  {
  }
  void ExampleProblem::evalObjGrad(RefVec)
  {
  }

  void ExampleProblem::evalLinCstr(RefVec, Index)
  {
  }
  void ExampleProblem::evalLinCstrGrad(RefVec, Index)
  {
  }

  void ExampleProblem::evalNonLinCstr(RefVec, Index)
  {
  }
  void ExampleProblem::evalNonLinCstrGrad(RefVec, Index)
  {
  }
}
