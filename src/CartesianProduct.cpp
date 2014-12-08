#include <pgsolver/CartesianProduct.h>

namespace pgs
{
  CartesianProduct::CartesianProduct()
    : Manifold(0,0)
  {
    startIndexT_.push_back(0);
    startIndexR_.push_back(0);
  }

  CartesianProduct::CartesianProduct(const Manifold& m1, const Manifold& m2)
    : Manifold(0,0)
  {
    startIndexT_.push_back(0);
    startIndexR_.push_back(0);
    multiply(m1);
    multiply(m2);
  }

  bool CartesianProduct::isValidInit_(const Eigen::VectorXd& val) const
  {
    bool out = true;
    for (std::size_t i = 0; i<numberOfSubmanifolds(); ++i)
    {
      out = out && submanifolds_[i]->isValidInit( 
          val.segment(startIndexR_[i], submanifolds_[i]->representationDim()));
    }
    return out;
  }

  CartesianProduct& CartesianProduct::multiply(const Manifold& m)
  {
    m.lock();
    setDimension(dim() + m.dim());
    setRepresentationDimension(representationDim() + m.representationDim());
    submanifolds_.push_back(&m);
    startIndexT_.push_back(startIndexT_.back() + m.dim());
    startIndexR_.push_back(startIndexR_.back() + m.representationDim());
    return *this;
  }

  size_t CartesianProduct::numberOfSubmanifolds() const
  {
    return submanifolds_.size();
  }

  const Manifold& CartesianProduct::operator()(size_t i) const
  {
    assert(i < submanifolds_.size() && "invalid index");
    return *submanifolds_[i];
  }

  Segment CartesianProduct::getValue(RefVec val, size_t i) const
  {
    assert(i < submanifolds_.size() && "invalid index");
    return val.segment(startIndexR_[i], submanifolds_[i]->representationDim());
  }

  ConstSegment CartesianProduct::getValueConst(const ConstRefVec& val, size_t i) const
  {
    assert(i < submanifolds_.size() && "invalid index");
    return val.segment(startIndexR_[i], submanifolds_[i]->representationDim());
  }

  std::string CartesianProduct::toString(const ConstRefVec& val, std::string& prefix) const
  {
    std::stringstream ss;
    std::string SubManPrefix("  ");
    for (std::size_t i = 0; i<numberOfSubmanifolds(); ++i)
    {
      ss << prefix << submanifolds_[i]->toString(
                          getValueConst(val, i), SubManPrefix);
    }
    return ss.str();
  }

  void CartesianProduct::plus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      submanifolds_[i]->plus(getValue(out, i), 
                              getValueConst(x, i), 
                              v.segment(startIndexT_[i], submanifolds_[i]->dim()));
    }
  }

  void CartesianProduct::minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      submanifolds_[i]->minus(out.segment(startIndexT_[i], 
                              submanifolds_[i]->dim()), 
                              getValueConst(x, i), 
                              getValueConst(y, i));
    }
  }

  void CartesianProduct::setIdentity_(RefVec out) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
      submanifolds_[i]->setIdentity(getValue(out, i));
  }

  Eigen::MatrixXd CartesianProduct::diffMap_(const ConstRefVec& x ) const
  {
    Eigen::MatrixXd J(representationDim(),dim());
    J.setZero();
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      J.block(startIndexR_[i],
              startIndexT_[i],
              submanifolds_[i]->representationDim(),
              submanifolds_[i]->dim()) 
        = submanifolds_[i]->diffMap(x.segment(startIndexR_[i], 
                                  submanifolds_[i]->representationDim()));
    }
    return J;
  }

  void CartesianProduct::applyDiffMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      Index startTi = startIndexT_[i];
      Index startRi = startIndexR_[i];
      Index repDimi = submanifolds_[i]->representationDim();
      Index dimi = submanifolds_[i]->dim();
      submanifolds_[i]->applyDiffMap( out.middleCols(startTi, dimi),
                                      in.middleCols(startRi, repDimi),
                                      x.segment(startRi, repDimi));
    }
  }
  
  Eigen::MatrixXd CartesianProduct::diffInvMap_(const ConstRefVec& x) const
  {
    Eigen::MatrixXd J(dim(),representationDim());
    J.setZero();
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      J.block(startIndexT_[i],
              startIndexR_[i],
              submanifolds_[i]->dim(),
              submanifolds_[i]->representationDim()) 
        = submanifolds_[i]->diffInvMap(x.segment(startIndexR_[i], 
                                  submanifolds_[i]->representationDim()));
    }
    return J;
  }

  void CartesianProduct::applyDiffInvMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      Index startTi = startIndexT_[i];
      Index startRi = startIndexR_[i];
      Index repDimi = submanifolds_[i]->representationDim();
      Index dimi = submanifolds_[i]->dim();
      submanifolds_[i]->applyDiffInvMap( out.middleCols(startRi, repDimi),
                                      in.middleCols(startTi, dimi),
                                      x.segment(startRi, repDimi));
    }
  }
}

