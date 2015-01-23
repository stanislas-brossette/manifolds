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
      static double LineSearch(Solver& s, 
                              Problem& p, ProblemEvaluation& probE,  
                              Filter& filter_, Eigen::VectorXd& step, 
                              SolverOptions& opti,
                              Eigen::VectorXi test = Eigen::VectorXi::Zero(0));
  };
}

#endif //_PGS_LINE_SEARCHER_H_
