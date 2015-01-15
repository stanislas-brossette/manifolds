namespace pgs
{
  enum eHessianUpdate { EXACT, BFGS, SR1 };
  struct OptimOptions
  {
    eHessianUpdate hessianUpdateMethod=BFGS; 
  };
}
