// Copyright (c) 2015 CNRS
// Authors: Stanislas Brossette, Adrien Escande

// This file is part of manifolds
// manifolds is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.

// manifolds is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// manifolds. If not, see
// <http://www.gnu.org/licenses/>.

#include <manifolds/CartesianProduct.h>
#include <manifolds/mnf_assert.h>

namespace mnf
{
  CartesianProduct::CartesianProduct()
    : Manifold(0,0,0)
  {
    startIndexT_.push_back(0);
    startIndexR_.push_back(0);
    name() = "";
  }

  CartesianProduct::CartesianProduct(std::initializer_list<Manifold*> m)
    : Manifold(0,0,0)
  {
    startIndexT_.push_back(0);
    startIndexR_.push_back(0);
    name() = "";
    for (auto mi : m)
      multiply(*(mi));
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

  void CartesianProduct::getIdentityOnTxM_(RefMat out, const ConstRefVec& x) const
  {
    for (std::size_t i = 0; i<numberOfSubmanifolds(); ++i)
    {
      submanifolds_[i]->getIdentityOnTxM(getView<T, T>(out,i), getConstView<R>(x,i));
      for (std::size_t j = 0; j<i; ++j)
      {
        getView<T, T>(out,i,j).setZero();
        getView<T, T>(out,j,i).setZero();
      }
    }
  }

  CartesianProduct& CartesianProduct::multiply(const Manifold& m)
  {
    testLock();
    m.lock();
    if(dim() != 0)
      name() += "x";
    name() += m.name();
    setDimension(dim() + m.dim());
    setTangentDimension(tangentDim() + m.tangentDim());
    setRepresentationDimension(representationDim() + m.representationDim());
    submanifolds_.push_back(std::shared_ptr<const Manifold>(copyManifold(m)));
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

  void CartesianProduct::display(std::string prefix) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      if(submanifolds_[i]->isElementary())
        std::cout << prefix << submanifolds_[i]->name()<< std::endl;
      else
      {
        std::cout << prefix << "/----------------------------------"<< std::endl;
        submanifolds_[i]->display(prefix + "| ");
        std::cout << prefix << "\\----------------------------------"<< std::endl;
      }
    }
  }

  const Manifold& CartesianProduct::operator()(size_t i) const
  {
    mnf_assert(i < submanifolds_.size() && "invalid index");
    return *submanifolds_[i];
  }

  std::string CartesianProduct::toString(const ConstRefVec& val, const std::string& prefix, const Eigen::IOFormat fmt) const
  {
    std::stringstream ss;
    size_t n = numberOfSubmanifolds();
    for (std::size_t i = 0; i<n-1; ++i)
    {
      ss << submanifolds_[i]->toString(getConstView<R>(val, i), prefix + "  ", fmt) << std::endl;;
    }
    ss << submanifolds_.back()->toString(getConstView<R>(val, n - 1), prefix + "  ", fmt);
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

  void CartesianProduct::applyInvTransportOnTheRight_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
    {
      submanifolds_[i]->applyInvTransportOnTheRight(getView<F, T>(out, i),
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

  bool CartesianProduct::isInTxM_(const ConstRefVec& x, const ConstRefVec& v, const double& prec) const
  {
    bool b = true;
    for (size_t i = 0; i < submanifolds_.size() && b; ++i)
      b = submanifolds_[i]->isInTxM(getConstView<R>(x, i), getConstView<T>(v, i), prec);
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

  void CartesianProduct::getTrustMagnitude_(RefVec out) const
  {
    for (size_t i = 0; i < submanifolds_.size(); ++i)
      submanifolds_[i]->getTrustMagnitude(getView<T>(out, i));
  }

  long CartesianProduct::getTypeId() const
  {
    return utils::hash::computeHash("CartesianProduct");
  }

  Manifold* CartesianProduct::getNewCopy() const
  {
    CartesianProduct* copy = new CartesianProduct(*this);
    copy->instanceId_ = this->instanceId_;

    return copy;
  }
}
