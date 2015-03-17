#include <manifolds/CartesianProduct.h>
#include <manifolds/pgs_assert.h>

namespace pgs
{
  CartesianProduct::CartesianProduct()
    : Manifold(0,0,0)
  {
    startIndexT_.push_back(0);
    startIndexR_.push_back(0);
    name() = "";
  }

  CartesianProduct::CartesianProduct(const Manifold& m1, const Manifold& m2)
    : Manifold(0,0,0)
  {
    startIndexT_.push_back(0);
    startIndexR_.push_back(0);
    multiply(m1);
    multiply(m2);
    name() = m1.name() + "x" + m2.name();
  }

  bool CartesianProduct::isInM_(const Eigen::VectorXd& val, const double& ) const
  {
    bool out = true;
    for (std::size_t i = 0; i<numberOfSubmanifolds(); ++i)
    {
      out = out && submanifolds_[i]->isInM(getConstView<R>(val, i));
    }
    return out;
  }
   
  void CartesianProduct::forceOnM_(RefVec out, const ConstRefVec& in) const
  {
    for (std::size_t i = 0; i<numberOfSubmanifolds(); ++i)
      submanifolds_[i]->forceOnM(getView<R>(out,i), getConstView<R>(in,i));
  }

  CartesianProduct& CartesianProduct::multiply(const Manifold& m)
  {
    m.lock();
    if(dim() != 0)
      name() += "x"; 
    name() += m.name();
    setDimension(dim() + m.dim());
    setTangentDimension(tangentDim() + m.tangentDim());
    setRepresentationDimension(representationDim() + m.representationDim());
    submanifolds_.push_back(&m);
    startIndexT_.push_back(startIndexT_.back() + m.tangentDim());
    startIndexR_.push_back(startIndexR_.back() + m.representationDim());
    return *this;
  }

  size_t CartesianProduct::numberOfSubmanifolds() const
  {
    return submanifolds_.size();
  }

  bool CartesianProduct::isElementary() const
  {
    return false;
  }

  const Manifold& CartesianProduct::operator()(size_t i) const
  {
    pgs_assert(i < submanifolds_.size() && "invalid index");
    return *submanifolds_[i];
  }

  std::string CartesianProduct::toString(const ConstRefVec& val, const std::string& prefix, int prec) const
  {
    std::stringstream ss;
    size_t n = numberOfSubmanifolds();
    for (std::size_t i = 0; i<n-1; ++i)
    {
      ss << submanifolds_[i]->toString(getConstView<R>(val, i), prefix + "  ", prec) << std::endl;;
    }
    ss << submanifolds_.back()->toString(getConstView<R>(val, n - 1), prefix + "  ", prec);
    return ss.str();
  }

  void CartesianProduct::createRandomPoint_(RefVec out, double coeff) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
      submanifolds_[i]->createRandomPoint(getView<R>(out, i), coeff);
  }

  void CartesianProduct::retractation_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      submanifolds_[i]->retractation(getView<R>(out, i),
                              getConstView<R>(x, i),
                              getConstView<T>(v, i));
    }
  }

  void CartesianProduct::pseudoLog_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      submanifolds_[i]->pseudoLog(getView<T>(out,i),
                              getConstView<R>(x, i),
                              getConstView<R>(y, i));
    }
  }

  void CartesianProduct::pseudoLog0_(RefVec out, const ConstRefVec& x) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      submanifolds_[i]->pseudoLog0(getView<T>(out,i),
                               getConstView<R>(x, i));
    }
  }

  void CartesianProduct::setZero_(RefVec out) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
      submanifolds_[i]->setZero(getView<R>(out, i));
  }

  Eigen::MatrixXd CartesianProduct::diffRetractation_(const ConstRefVec& x ) const
  {
    Eigen::MatrixXd J(representationDim(),tangentDim());
    J.setZero();
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      getView<R, T>(J, i) = submanifolds_[i]->diffRetractation(getConstView<R>(x, i));
    }
    return J;
  }

  void CartesianProduct::applyDiffRetractation_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      submanifolds_[i]->applyDiffRetractation( getView<F, T>(out, i),
                                      getConstView<F, R>(in, i),
                                      getConstView<R>(x, i));
    }
  }

  Eigen::MatrixXd CartesianProduct::diffPseudoLog0_(const ConstRefVec& x) const
  {
    Eigen::MatrixXd J(tangentDim(),representationDim());
    J.setZero();
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      getView<T, R>(J, i) = submanifolds_[i]->diffPseudoLog0(getConstView<R>(x, i));
    }
    return J;
  }

  void CartesianProduct::applyDiffPseudoLog0_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      submanifolds_[i]->applyDiffPseudoLog0(getView<F,R>(out, i),
                                        getConstView<F, T>(in, i),
                                        getConstView<R>(x, i));
    }
  }

  void CartesianProduct::applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      submanifolds_[i]->applyTransport(getView<T, F>(out, i),
                                       getConstView<T, F>(in, i),
                                       getConstView<R>(x, i),
                                       getConstView<T>(v, i));
    }
  }

  void CartesianProduct::applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      submanifolds_[i]->applyInvTransport(getView<F, T>(out, i),
                                          getConstView<F, T>(in, i),
                                          getConstView<R>(x, i),
                                          getConstView<T>(v, i));
    }
  }

  void CartesianProduct::tangentConstraint_(RefMat out, const ConstRefVec& x) const
  {
    Index k = 0;
    out.setZero();
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      Index s = submanifolds_[i]->tangentDim() - submanifolds_[i]->dim();
      auto Ci = getView<F, T>(out, i);
      submanifolds_[i]->tangentConstraint(Ci.middleRows(k,s),
                                          getConstView<R>(x, i));
      k += s;
    }
  }

  bool CartesianProduct::isInTxM_(const ConstRefVec& x, const ConstRefVec& v, const double& ) const
  {
    bool b = true;
    for (size_t i = 0; i < submanifolds_.size() && b; ++i)
      b = submanifolds_[i]->isInTxM(getConstView<R>(x, i), getConstView<T>(v, i));
    return b;
  }

  void CartesianProduct::forceOnTxM_(RefVec out, const ConstRefVec& in, const ConstRefVec& x) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
      submanifolds_[i]->forceOnTxM(getView<T>(out, i), 
                                   getConstView<T>(in, i),
                                   getConstView<R>(x, i));
  }

  void CartesianProduct::limitMap_(RefVec out) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
      submanifolds_[i]->limitMap(getView<T>(out, i)); 
  }

  void CartesianProduct::getTypicalMagnitude_(RefVec out) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
      submanifolds_[i]->getTypicalMagnitude(getView<T>(out, i)); 
  }
}

