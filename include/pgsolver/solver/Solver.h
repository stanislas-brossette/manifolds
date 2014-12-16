#ifndef _PGS_SOLVER_H_
#define _PGS_SOLVER_H_

#include <iostream>
#include <Eigen/Core>

#include <pgsolver/solver/Problem.h>
#include <pgsolver/solver/Results.h>
#include <pgsolver/manifolds/Point.h>
#include <pgsolver/solver/TempData.h>
#include <pgsolver/solver/ConstraintManager.h>

namespace pgs
{
  class Solver
  {
    public:
      Solver();
      Results solve(Problem& problem, Point& x0);

    private:
      void initSolver(Problem& problem, Point& x0);

    private:
      ProblemEvaluation probEval_;
      SolverOption opt_;
      ConstraintManager cstrMngr_;
  };
}

#endif //_PGS_SOLVER_H_
