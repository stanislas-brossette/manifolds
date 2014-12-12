#ifndef _PGS_PROBLEM_ON_REPRESENTATION_SPACE_H_
#define _PGS_PROBLEM_ON_REPRESENTATION_SPACE_H_

#include <iostream>
#include <Eigen/Core>

#include <pgsolver/manifolds/defs.h>
#include <pgsolver/solver/Problem.h>

namespace pgs
{
  /// \brief This kind of problem has its objective and constraints written on
  /// the representation space of the manifold and does the compositions with the
  /// map automatically
  class ProblemOnRepresentationSpace : public Problem
  {
    public:
      ProblemOnRepresentationSpace(Manifold& M);
      ProblemOnRepresentationSpace(Manifold& manifold, const Point& x);
      virtual void getUB(RefVec out) const;
      virtual void getLB(RefVec out) const;

      virtual void getCstrLB(RefVec out, size_t i) const;
      virtual void getCstrUB(RefVec out, size_t i) const;

      virtual void evalObj(double& out) const;
      virtual void evalObjGrad(RefMat out) const;

      virtual void evalLinCstr(RefVec out, size_t i) const;
      virtual void evalLinCstrGrad(RefVec out, size_t i) const;

      virtual void evalNonLinCstr(RefVec out, size_t i) const;
      virtual void evalNonLinCstrGrad(RefVec out, size_t i) const;
  };
}

#endif //_PGS_PROBLEM_ON_REPRESENTATION_SPACE_H_
