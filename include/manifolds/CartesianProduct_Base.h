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

#ifndef _MANIFOLDS_CARTESIAN_PRODUCT_H_
#define _MANIFOLDS_CARTESIAN_PRODUCT_H_

#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <memory>

#include <manifolds/defs.h>
#include <manifolds/Manifold_Base.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/utils.h>

namespace mnf
{
  /// \brief Manifold representing the cartesian product of several submanifolds
  class CartesianProduct_Base : public Manifold_Base
  {
    friend CartesianProduct;
  public:
    /// \brief Default constructor
    CartesianProduct_Base();

    /// \brief Constructor of the manifold from a list of manifolds
    CartesianProduct_Base(const std::initializer_list<Manifold*> m);

    /// \brief Constructor of the manifold composed of \f$ m1\times m2\f$
    CartesianProduct_Base(const Manifold& m1, const Manifold& m2);

    /// \brief Adds manifold m to the current composed manifold\n
    /// This method cannot be executed if the manifold is locked
    CartesianProduct_Base& multiply(const Manifold& m);

    virtual size_t numberOfSubmanifolds() const;
    virtual const Manifold& operator()(const size_t i) const;

    virtual std::string toString(const ConstRefVec& val, const std::string& prefix = "", const Eigen::IOFormat& fmt = mnf::defaultFormat) const;

    virtual bool isElementary() const;

    virtual void display(const std::string& prefix = "") const;

  protected:
    virtual bool isInM_(const Eigen::VectorXd& val, double prec) const;
    virtual void forceOnM_(RefVec out, const ConstRefVec& in) const;
    virtual void getIdentityOnTxM_(RefMat out, const ConstRefVec& x) const;

    virtual Index startR(const size_t i) const;
    virtual Index startT(const size_t i) const;

    virtual void createRandomPoint_(RefVec out, const double coeff) const;
    virtual void retractation_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void pseudoLog_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const;
    virtual void pseudoLog0_(RefVec out, const ConstRefVec& x) const;
    virtual void setZero_(RefVec out) const;
    virtual Eigen::MatrixXd diffRetractation_(const ConstRefVec& x) const;
    virtual void applyDiffRetractation_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    virtual Eigen::MatrixXd diffPseudoLog0_(const ConstRefVec& x) const;
    virtual void applyDiffPseudoLog0_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    virtual void applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void applyInvTransportOnTheRight_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const;

    virtual void tangentConstraint_(RefMat out, const ConstRefVec& x) const;
    virtual bool isInTxM_(const ConstRefVec& x, const ConstRefVec& v, const double& prec) const;
    virtual void forceOnTxM_(RefVec out, const ConstRefVec& in, const ConstRefVec& x) const;
    virtual void limitMap_(RefVec out) const;
    virtual void getTypicalMagnitude_(RefVec out) const;
    virtual void getTrustMagnitude_(RefVec out) const;
    virtual long getTypeId() const;

    virtual Manifold_ptr getNewCopy() const;

  private:
    /// \brief List of pointers on all the manifolds in the cartesian product
    std::vector<std::shared_ptr<const Manifold>> submanifolds_;

    /// \brief List of start index of submanifolds in a vector of the
    /// tangent space
    std::vector<Index> startIndexT_;

    /// \brief List of start index of submanifolds in a vector of the
    /// representation space
    std::vector<Index> startIndexR_;
  };

  inline Index CartesianProduct_Base::startR(const size_t i) const
  {
    mnf_assert(i < numberOfSubmanifolds() && "invalid index");
    return startIndexR_[i];
  }

  inline Index CartesianProduct_Base::startT(const size_t i) const
  {
    mnf_assert(i < numberOfSubmanifolds() && "invalid index");
    return startIndexT_[i];
  }
}

#endif //_MANIFOLDS_CARTESIAN_PRODUCT_H_
