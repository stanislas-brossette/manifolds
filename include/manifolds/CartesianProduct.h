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
#include <manifolds/defs.h>
#include <manifolds/Manifold.h>
#include <manifolds/utils.h>

namespace mnf
{
  /// \brief Manifold representing the cartesian product of several submanifolds
  class MANIFOLDS_API CartesianProduct : public Manifold
  {
  public:
    /// \brief Default constructor
    CartesianProduct();

    /// \brief Constructor of the manifold from a list of manifolds
    CartesianProduct(std::initializer_list<Manifold*> m);

    /// \brief Constructor of the manifold composed of \f$ m1\times m2\f$
    CartesianProduct(const Manifold& m1, const Manifold& m2);

    /// \brief Adds manifold m to the current composed manifold\n
    /// This method cannot be executed if the manifold is locked
    CartesianProduct& multiply(const Manifold& m);

    virtual size_t numberOfSubmanifolds() const;
    virtual const Manifold& operator()(size_t i) const;

    virtual std::string toString(const ConstRefVec& val, const std::string& prefix = "", int prec = 6) const;

    virtual bool isElementary() const;

    virtual void display(std::string prefix = "") const;

  protected:
    virtual bool isInM_(const Eigen::VectorXd& val, const double& prec) const;
    virtual void forceOnM_(RefVec out, const ConstRefVec& in) const;

    virtual Index startR(size_t i) const;
    virtual Index startT(size_t i) const;

    virtual void createRandomPoint_(RefVec out, double coeff) const;
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
    virtual long getTypeId() const;

  private:
    /// \brief List of pointers on all the manifolds in the cartesian product
    std::vector<const Manifold* > submanifolds_;

    /// \brief List of start index of submanifolds in a vector of the
    /// tangent space
    std::vector<Index> startIndexT_;

    /// \brief List of start index of submanifolds in a vector of the
    /// representation space
    std::vector<Index> startIndexR_;
  };

  inline Index CartesianProduct::startR(size_t i) const
  {
    mnf_assert(i < numberOfSubmanifolds() && "invalid index");
    return startIndexR_[i];
  }

  inline Index CartesianProduct::startT(size_t i) const
  {
    mnf_assert(i < numberOfSubmanifolds() && "invalid index");
    return startIndexT_[i];
  }
}

#endif //_MANIFOLDS_CARTESIAN_PRODUCT_H_
