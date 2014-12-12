#include <pgsolver/solver/ProblemOnRepresentationSpace.h>

namespace pgs
{
      
  ProblemOnRepresentationSpace::ProblemOnRepresentationSpace(Manifold& M)
    : Problem(M)
  {
  }
  ProblemOnRepresentationSpace::ProblemOnRepresentationSpace(Manifold& M, const Point& x)
    : Problem(M, x)
  {
  }

  void ProblemOnRepresentationSpace::getUB(RefVec) const
  {
  }
  void ProblemOnRepresentationSpace::getLB(RefVec) const
  {
  }

  void ProblemOnRepresentationSpace::getCstrLB(RefVec, size_t) const
  {
  }
  void ProblemOnRepresentationSpace::getCstrUB(RefVec, size_t) const
  {
  }

  void ProblemOnRepresentationSpace::evalObj(double&) const
  {
  }
  void ProblemOnRepresentationSpace::evalObjGrad(RefMat) const
  {
  }

  void ProblemOnRepresentationSpace::evalLinCstr(RefVec, size_t) const
  {
  }
  void ProblemOnRepresentationSpace::evalLinCstrGrad(RefVec, size_t) const
  {
  }

  void ProblemOnRepresentationSpace::evalNonLinCstr(RefVec, size_t) const
  {
  }
  void ProblemOnRepresentationSpace::evalNonLinCstrGrad(RefVec, size_t) const
  {
  }
}
