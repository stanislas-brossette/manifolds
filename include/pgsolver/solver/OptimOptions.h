#ifndef _PGS_OPTIM_OPTIONS_H_
#define _PGS_OPTIM_OPTIONS_H_

namespace pgs
{
  enum eHessianUpdate { EXACT, BFGS, SR1 };
  struct OptimOptions
  {
    eHessianUpdate hessianUpdateMethod=BFGS; 
  };
}
#endif //_PGS_OPTIM_OPTIONS_H_
