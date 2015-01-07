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
      virtual void evalObjDiff(RefMat out) const;

      virtual void evalLinCstr(RefMat out, size_t i) const;
      virtual void evalLinCstrDiff(RefMat out, size_t i) const;
      virtual void getLinCstrLB(RefMat out, size_t i) const;
      virtual void getLinCstrUB(RefMat out, size_t i) const;

      virtual void evalNonLinCstr(RefMat out, size_t i) const;
      virtual void evalNonLinCstrDiff(RefMat out, size_t i) const;
      virtual void getNonLinCstrLB(RefMat out, size_t i) const;
      virtual void getNonLinCstrUB(RefMat out, size_t i) const;

      virtual size_t numberOfCstr() const;
      virtual Index linCstrDim(size_t i) const;
      virtual Index nonLinCstrDim(size_t i) const;

    private:
      static RealSpace R3;

      double a, b, c, d; //Coefficient of the plan for the linear cstr
      double R1, R2; //Inner and outer radiusesof the spheric envelope
  };
}

#endif //_PGS_EXAMPLE_GEOMETRIC_PROBLEM_H_
