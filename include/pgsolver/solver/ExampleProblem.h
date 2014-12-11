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
      ExampleProblem(const Manifold& M);
      virtual void getUB(RefVec out);
      virtual void getLB(RefVec out);

      virtual void getCstrLB(RefVec out);
      virtual void getCstrUB(RefVec out);

      virtual void evalObj(RefVec out);
      virtual void evalObjGrad(RefVec out);

      virtual void evalLinCstr(RefVec out, Index i);
      virtual void evalLinCstrGrad(RefVec out, Index i);

      virtual void evalNonLinCstr(RefVec out, Index i);
      virtual void evalNonLinCstrGrad(RefVec out, Index i);
  };
}

#endif //_PGS_EXAMPLE_PROBLEM_H_
