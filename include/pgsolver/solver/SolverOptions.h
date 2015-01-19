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
  struct SolverOptions
  {
    int maxIter = 100;
    double epsilon_P = 1e-3;
    double epsilon_D = 1e-3;
    double gammaFilter = 1e-16;
    Filter::eOption filterOpt = Filter::EXISTING;

    //TODO: make a pair of {hessMeth, hessType}
    eHessianUpdateMethod hessianUpdateMethod = BFGS; 
    eHessianUpdateType hessianUpdateType = GROUPED; 
    eGlobalization globalizationMethod = LINESEARCH;
    eLineSearchType lineSearchMethod = FILTER;
  };
}
#endif //_PGS_SOLVER_OPTIONS_H_
