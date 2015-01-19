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
    totalDimLin_ = 0;
    totalDimNonLin_ = 0;
    for (size_t i = 0; i<nbCstr; ++i)
    {
      dimLin_[i] = p.linCstrDim(i);
      dimNonLin_[i] = p.nonLinCstrDim(i);
      startLin_[i] = totalDimLin_;
      startNonLin_[i] = totalDimNonLin_;
      totalDimLin_ += dimLin_[i];
      totalDimNonLin_ += dimNonLin_[i];
    }
    totalDim_ = totalDimLin_ + totalDimNonLin_;
  }

  RefMat ConstraintManager::getViewLin(RefMat J, size_t i) const
  {
    assert(J.rows() == totalDimLin_ && "Wrong nunmber of lines");
    return J.middleRows(startLin_[i],dimLin_[i]);
  }
  const ConstRefMat ConstraintManager::getConstViewLin(const ConstRefMat J, size_t i) const
  {
    assert(J.rows() == totalDimLin_ && "Wrong nunmber of lines");
    return J.middleRows(startLin_[i],dimLin_[i]);
  }

  RefMat ConstraintManager::getViewNonLin(RefMat J, size_t i) const
  {
    assert(J.rows() == totalDimNonLin_ && "Wrong nunmber of lines");
    return J.middleRows(startNonLin_[i],dimNonLin_[i]);
  }
  const ConstRefMat ConstraintManager::getConstViewNonLin(const ConstRefMat J, size_t i) const
  {
    assert(J.rows() == totalDimNonLin_ && "Wrong nunmber of lines");
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
  const Index& ConstraintManager::totalDim() const
  {
    return totalDim_;
  }
}
