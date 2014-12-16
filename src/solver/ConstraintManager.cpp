#include <Eigen/Core>
#include <pgsolver/manifolds/defs.h>
#include <pgsolver/solver/ConstraintManager.h>

namespace pgs
{
  ConstraintManager::ConstraintManager()
  {
    std::cout << "Constraint Manager"<< std::endl;
  }
  
  void ConstraintManager::init(Problem& p)
  {
    size_t nbCstr = p.numberOfCstr();
    std::cout << "nbCstr = " << nbCstr << std::endl;
    startLin.resize(nbCstr);
    dimLin.resize(nbCstr);
    startNonLin.resize(nbCstr);
    dimNonLin.resize(nbCstr);

    Index startIndex = 0;
    for (size_t i = 0; i<nbCstr; ++i)
    {
      dimLin[i] = p.linCstrDim(i);
      dimNonLin[i] = p.nonLinCstrDim(i);
      startLin[i] = startIndex;
      startIndex += dimLin[i];
    }
    for (size_t i = 0; i<nbCstr; ++i)
    {
      startNonLin[i] = startIndex;
      startIndex += dimNonLin[i];
    }
  }
  RefMat ConstraintManager::getViewLin(RefMat J, size_t i)
  {
    return J.middleRows(startLin[i],dimLin[i]);
  }
  RefMat ConstraintManager::getViewNonLin(RefMat J, size_t i)
  {
    return J.middleRows(startNonLin[i],dimNonLin[i]);
  }
}
