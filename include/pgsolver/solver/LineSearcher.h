#ifndef _PGS_LINE_SEARCHER_H_
#define _PGS_LINE_SEARCHER_H_

#include <Eigen/Core>

#include <pgsolver/solver/Problem.h>
#include <pgsolver/solver/Solver.h>
#include <pgsolver/solver/SolverOptions.h>

namespace pgs
{
  class LineSearcher
  {
    public:
      static double LineSearch(Solver& s, Problem& p, Filter& filter_, Eigen::VectorXd& step);
  };
}

#endif //_PGS_LINE_SEARCHER_H_
