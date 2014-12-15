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

      virtual void evalObj(double& out) const;
      virtual void evalObjGrad(RefMat out) const;

      virtual void evalLinCstr(RefVec out, size_t i) const;
      virtual void evalLinCstrGrad(RefMat out, size_t i) const;
      virtual void getLinCstrLB(RefVec out, size_t i) const;
      virtual void getLinCstrUB(RefMat out, size_t i) const;

      virtual void evalNonLinCstr(RefVec out, size_t i) const;
      virtual void evalNonLinCstrGrad(RefMat out, size_t i) const;
      virtual void getNonLinCstrLB(RefVec out, size_t i) const;
      virtual void getNonLinCstrUB(RefVec out, size_t i) const;
    private:
      static RealSpace R3;

      double a, b, c, d; //Coefficient of the plan for the linear cstr
      double R1, R2; //Inner and outer radiusesof the spheric envelope

  };
}

#endif //_PGS_EXAMPLE_GEOMETRIC_PROBLEM_H_
