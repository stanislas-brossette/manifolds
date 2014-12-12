#ifndef _PGS_EXAMPLE_PROBLEM_H_
#define _PGS_EXAMPLE_PROBLEM_H_

#include <iostream>
#include <Eigen/Core>

#include <pgsolver/manifolds/defs.h>
#include <pgsolver/solver/Problem.h>

namespace pgs
{
  class ExampleProblem : public Problem
  {
    public:
      ExampleProblem(Manifold& M);
      ExampleProblem(Manifold& M, const Point& x);
      virtual void getUB(RefVec out) const;
      virtual void getLB(RefVec out) const;

      virtual void getCstrLB(RefVec out, size_t i) const;
      virtual void getCstrUB(RefVec out, size_t i) const;

      virtual void evalObj(RefVec out) const;
      virtual void evalObjGrad(RefVec out) const;

      virtual void evalLinCstr(RefVec out, size_t i) const;
      virtual void evalLinCstrGrad(RefVec out, size_t i) const;

      virtual void evalNonLinCstr(RefVec out, size_t i) const;
      virtual void evalNonLinCstrGrad(RefVec out, size_t i) const;
  };
}

#endif //_PGS_EXAMPLE_PROBLEM_H_
