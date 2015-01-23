#ifndef _PGS_SOLVER_OPTIONS_H_
#define _PGS_SOLVER_OPTIONS_H_

#include <pgsolver/solver/Filter.h>

namespace pgs
{
  enum eHessianUpdateMethod { EXACT, BFGS, SR1 };
  enum eHessianUpdateType { GROUPED, INDIVIDUAL };
  enum eGlobalization { NONE, LINESEARCH, TRUSTREGION };
  enum eLineSearchType { FILTER };
  enum eTrustRegionType { };
  enum eCstrStatus {
    VIOLATED_LB = -2,
    VIOLATED_UB = -1,
    SATISFIED = 0,
    ACTIVE_LB = 1,
    ACTIVE_UB = 2,
    ACTIVE_EQUALITY = 3
    };
  struct SolverOptions
  {
    int maxIter = 100;
    double epsilon_P = 1e-3;
    double epsilon_D = 1e-3;
    double gammaFilter = 1e-16;
    double epsilonFeasibility = 1e-3;
    Filter::eOption filterOpt = Filter::EXISTING;
    int VERBOSE = 1;

    //TODO: make a pair of {hessMeth, hessType}
    eHessianUpdateMethod hessianUpdateMethod = BFGS;
    eHessianUpdateType hessianUpdateType = GROUPED;
    eGlobalization globalizationMethod = LINESEARCH;
    eLineSearchType lineSearchMethod = FILTER;
  };
}
#endif //_PGS_SOLVER_OPTIONS_H_
