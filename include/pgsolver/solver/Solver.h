#ifndef _PGS_SOLVER_H_
#define _PGS_SOLVER_H_

#include <iostream>
#include <Eigen/Core>

#include <pgsolver/solver/Problem.h>
#include <pgsolver/manifolds/Point.h>

namespace pgs
{
  class Solver
  {
    public:
      void solve(Problem& problem, Point& x0);
  };
}

#endif //_PGS_SOLVER_H_
