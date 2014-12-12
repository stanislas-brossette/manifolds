#ifndef _PGS_EXAMPLE_GEOMETRIC_PROBLEM_H_
#define _PGS_EXAMPLE_GEOMETRIC_PROBLEM_H_

#include <iostream>
#include <Eigen/Core>

#include <pgsolver/manifolds/defs.h>
#include <pgsolver/solver/Problem.h>
#include <pgsolver/manifolds/RealSpace.h>

namespace pgs
{
  class ExampleGeometricProblem : public Problem
  {
    public:
      ExampleGeometricProblem();

      virtual void getTangentLB(RefVec out) const;
      virtual void getTangentUB(RefVec out) const;

      virtual void evalObj(RefVec out) const;
      virtual void evalObjGrad(RefVec out) const;

      virtual void evalLinCstr(RefVec out, size_t i) const;
      virtual void evalLinCstrGrad(RefVec out, size_t i) const;
      virtual void getLinCstrLB(RefVec out, size_t i) const;
      virtual void getLinCstrUB(RefVec out, size_t i) const;

      virtual void evalNonLinCstr(RefVec out, size_t i) const;
      virtual void evalNonLinCstrGrad(RefVec out, size_t i) const;
      virtual void getNonLinCstrLB(RefVec out, size_t i) const;
      virtual void getNonLinCstrUB(RefVec out, size_t i) const;
  };
}

#endif //_PGS_EXAMPLE_GEOMETRIC_PROBLEM_H_
