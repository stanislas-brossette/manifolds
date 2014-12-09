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
      out = out && submanifolds_[i]->isValidInit(getValueConst(val, i)); 
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

  //Segment CartesianProduct::getValue(RefVec val, size_t i) const
  //{
  //  assert(val.size() == representationDim());
  //  assert(i < submanifolds_.size() && "invalid index");
  //  return val.segment(startIndexR_[i], submanifolds_[i]->representationDim());
  //}

  //ConstSegment CartesianProduct::getValueConst(const ConstRefVec& val, size_t i) const
  //{
  //  assert(val.size() == representationDim());
  //  assert(i < submanifolds_.size() && "invalid index");
  //  return val.segment(startIndexR_[i], submanifolds_[i]->representationDim());
  //}

  //Segment CartesianProduct::getValueTangent(RefVec val, size_t i) const
  //{
  //  assert(val.size() == dim());
  //  assert(i < submanifolds_.size() && "invalid index");
  //  return val.segment(startIndexT_[i], submanifolds_[i]->dim());
  //}

  //ConstSegment CartesianProduct::getValueTangentConst(const ConstRefVec& val, size_t i) const
  //{
  //  assert(val.size() == dim());
  //  assert(i < submanifolds_.size() && "invalid index");
  //  return val.segment(startIndexT_[i], submanifolds_[i]->dim());
  //}

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
                              getValueTangentConst(v, i));
    }
  }

  void CartesianProduct::minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      submanifolds_[i]->minus(getValueTangent(out,i),
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
      //J.block(startIndexR_[i],
      //        startIndexT_[i],
      //        submanifolds_[i]->representationDim(),
      //        submanifolds_[i]->dim()) 
      //  = submanifolds_[i]->diffMap(getValueConst(x, i));
      getView<R, T>(J, i) = submanifolds_[i]->diffMap(getValueConst(x, i));
    }
    return J;
  }

  void CartesianProduct::applyDiffMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      submanifolds_[i]->applyDiffMap( getView<F, T>(out, i),
                                      getConstView<F, R>(in,i),
                                      getValueConst(x,i));
    }
  }
  
  Eigen::MatrixXd CartesianProduct::diffInvMap_(const ConstRefVec& x) const
  {
    Eigen::MatrixXd J(dim(),representationDim());
    J.setZero();
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      //J.block(startIndexT_[i],
      //        startIndexR_[i],
      //        submanifolds_[i]->dim(),
      //        submanifolds_[i]->representationDim()) 
      //  = submanifolds_[i]->diffInvMap(getValueConst(x,i));
      getView<T, R>(J, i) = submanifolds_[i]->diffInvMap(getValueConst(x, i));
    }
    return J;
  }

  void CartesianProduct::applyDiffInvMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      submanifolds_[i]->applyDiffInvMap(getView<F,R>(out, i),
                                        getConstView<F, T>(in,i),
                                        getValueConst(x,i));
    }
  }
}

