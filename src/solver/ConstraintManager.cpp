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
    startLin_.resize(nbCstr);
    dimLin_.resize(nbCstr);
    startNonLin_.resize(nbCstr);
    dimNonLin_.resize(nbCstr);
    Index startIndex = 0;
    for (size_t i = 0; i<nbCstr; ++i)
    {
      dimLin_[i] = p.linCstrDim(i);
      dimNonLin_[i] = p.nonLinCstrDim(i);
      startLin_[i] = startIndex;
      startIndex += dimLin_[i];
    }
    totalDimLin_ = startIndex;
    for (size_t i = 0; i<nbCstr; ++i)
    {
      startNonLin_[i] = startIndex;
      startIndex += dimNonLin_[i];
    }
    totalDimNonLin_ = startIndex - totalDimLin_;
  }

  RefMat ConstraintManager::getViewLin(RefMat J, size_t i)
  {
    return J.middleRows(startLin_[i],dimLin_[i]);
  }
  const ConstRefMat ConstraintManager::getViewLin(const ConstRefMat J, size_t i) const
  {
    return J.middleRows(startLin_[i],dimLin_[i]);
  }

  RefMat ConstraintManager::getViewNonLin(RefMat J, size_t i)
  {
    return J.middleRows(startNonLin_[i],dimNonLin_[i]);
  }
  const ConstRefMat ConstraintManager::getViewNonLin(const ConstRefMat J, size_t i) const
  {
    return J.middleRows(startNonLin_[i],dimNonLin_[i]);
  }

  const Index& ConstraintManager::totalDimLin() const
  {
    return totalDimLin_;
  }
  const Index& ConstraintManager::totalDimNonLin() const
  {
    return totalDimNonLin_;
  }
}
