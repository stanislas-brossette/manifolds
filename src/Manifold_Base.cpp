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

#include <stdexcept>
#include <manifolds/Manifold_Base.h>
#include <manifolds/Point.h>
#include <manifolds/mnf_assert.h>

namespace mnf
{
  long Manifold_Base::manifoldCounter_ = 0;

  Manifold_Base::Manifold_Base(Index dimension, Index tangentDimension, Index representationDimension)
    : dimension_(dimension)
    , tangentDim_(tangentDimension)
    , representationDim_(representationDimension)
    , lock_(false)
  {
    mnf_assert(0 <= dimension && "Negative dimension not accepted");
    mnf_assert(dimension <= tangentDimension);
    mnf_assert(tangentDimension <= representationDimension);
    trustMagnitude_.resize(tangentDimension);
    typicalMagnitude_.resize(tangentDimension);
    this->instanceId_ = manifoldCounter_++;
  }

  Manifold_Base::~Manifold_Base()
  {
  }

  Point Manifold_Base::createPoint() const
  {
    lock();
    return Point(ConstManifold(shared_from_this()));
  }

  Point Manifold_Base::createPoint(const ConstRefVec& val) const
  {
    if (isInM(val))
    {
      lock();
      return Point(shared_from_this(), val);
    }
    else
    {
      throw std::runtime_error("Bad Point Initialization");
    }
  }

  bool Manifold_Base::isInM(const Eigen::VectorXd& val, double prec) const
  {
    mnf_assert(val.size() == representationDim());
    return isInM_(val, prec);
  }

  void Manifold_Base::forceOnM(RefVec out, const ConstRefVec& in) const
  {
    mnf_assert(out.size() == representationDim());
    mnf_assert(in.size() == representationDim());
    return forceOnM_(out, in);
  }


  Point Manifold_Base::getZero() const
  {
    lock();
    Eigen::VectorXd id(representationDim_);
    setZero(id);
    return Point(shared_from_this(), id);
  }

  Point Manifold_Base::createRandomPoint(double coeff) const
  {
    lock();
    Eigen::VectorXd val(representationDim_);
    createRandomPoint(val, coeff);
    return Point(shared_from_this(), val);
  }

  void Manifold_Base::createRandomPoint(RefVec out, double coeff) const
  {
    mnf_assert(out.size() == representationDim_ && "wrong dimension in Manifold_Base::createRandomPoint");
    createRandomPoint_(out, coeff);
  }

  Index Manifold_Base::dim() const
  {
    return dimension_;
  }

  Index Manifold_Base::tangentDim() const
  {
    return tangentDim_;
  }

  Index Manifold_Base::representationDim() const
  {
    return representationDim_;
  }

  void Manifold_Base::display(const std::string& prefix) const
  {
    std::cout << prefix << name() << std::endl;
  }

  void Manifold_Base::retractation(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    mnf_assert(out.size() == representationDim_);
    mnf_assert(x.size() == representationDim_);
    mnf_assert(v.size() == tangentDim_);
    mnf_assert(isInTxM(x, v) && "Wrong tangent vector provided to retractation");
    retractation_(out, x, v);
  }

  void Manifold_Base::pseudoLog(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    mnf_assert(out.size() == tangentDim_);
    mnf_assert(x.size() == representationDim_);
    mnf_assert(y.size() == representationDim_);
    pseudoLog_(out, x, y);
  }

  void Manifold_Base::pseudoLog0(RefVec out, const ConstRefVec& x) const
  {
    mnf_assert(out.size() == tangentDim_);
    mnf_assert(x.size() == representationDim_);
    pseudoLog0_(out, x);
  }

  Eigen::MatrixXd Manifold_Base::diffRetractation(const ConstRefVec& x) const
  {
    mnf_assert(x.size() == representationDim_);
    return diffRetractation_(x);
  }

  void Manifold_Base::applyDiffRetractation(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    mnf_assert(in.cols() == representationDim_);
    mnf_assert(out.cols() == tangentDim_);
    mnf_assert(in.rows() == out.rows());
    mnf_assert(x.size() == representationDim_);
    applyDiffRetractation_(out, in, x);
  }

  Eigen::MatrixXd Manifold_Base::diffPseudoLog0(const ConstRefVec& x) const
  {
    mnf_assert(x.size() == representationDim_);
    return diffPseudoLog0_(x);
  }

  void Manifold_Base::applyDiffPseudoLog0(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    mnf_assert(out.cols() == representationDim_);
    mnf_assert(in.cols() == tangentDim_);
    mnf_assert(in.rows() == out.rows());
    mnf_assert(x.size() == representationDim_);
    applyDiffPseudoLog0_(out, in, x);
  }

  void Manifold_Base::applyTransport(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    mnf_assert(in.rows() == tangentDim_);
    mnf_assert(out.rows() == tangentDim_);
    mnf_assert(in.cols() == out.cols());
    mnf_assert(x.size() == representationDim());
    mnf_assert(v.size() == tangentDim_);
    mnf_assert(isInTxM(x, v));

    applyTransport_(out, in, x, v);
  }

  void Manifold_Base::applyInvTransport(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    mnf_assert(in.rows() == tangentDim_);
    mnf_assert(out.rows() == tangentDim_);
    mnf_assert(in.cols() == out.cols());
    mnf_assert(x.size() == representationDim());
    mnf_assert(v.size() == tangentDim_);
    mnf_assert(isInTxM(x, v));

    applyInvTransport_(out, in, x, v);
  }

  void Manifold_Base::applyInvTransportOnTheRight(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    mnf_assert(in.cols() == tangentDim_);
    mnf_assert(out.cols() == tangentDim_);
    mnf_assert(in.rows() == out.rows());
    mnf_assert(x.size() == representationDim());
    mnf_assert(v.size() == tangentDim_);
    mnf_assert(isInTxM(x, v));

    applyInvTransportOnTheRight_(out, in, x, v);
  }

  void Manifold_Base::setDimension(Index d)
  {
    mnf_assert(d>=0 && "Negative dimension not accepted");
    testLock();
    dimension_ = d;
  }

  void Manifold_Base::setTangentDimension(Index td)
  {
    mnf_assert(td >= 0 && "Negative dimension not accepted");
    testLock();
    tangentDim_ = td;
  }

  void Manifold_Base::setRepresentationDimension(Index rd)
  {
    mnf_assert(rd>=0 && "Negative dimension not accepted");
    testLock();
    representationDim_ = rd;
  }


  void Manifold_Base::setZero(RefVec out) const
  {
    mnf_assert(out.size() == representationDim_);
    setZero_(out);
  }

  void Manifold_Base::tangentConstraint(RefMat out, const ConstRefVec& x) const
  {
    mnf_assert(out.rows() == tangentDim_ - dimension_);
    mnf_assert(out.cols() == tangentDim_);
    mnf_assert(x.size() == representationDim());
    tangentConstraint_(out, x);
  }

  bool Manifold_Base::isInTxM(const ConstRefVec& x, const ConstRefVec& v, const double& prec) const
  {
    mnf_assert(v.size() == tangentDim_);
    mnf_assert(x.size() == representationDim());
    return isInTxM_(x, v, prec);
  }

  void Manifold_Base::forceOnTxM(RefVec out, const ConstRefVec& in, const ConstRefVec&x) const
  {
    mnf_assert(out.size() == tangentDim_);
    mnf_assert(x.size() == representationDim());
    mnf_assert(in.size() == tangentDim_);
    forceOnTxM_(out, in, x);
  }

  void Manifold_Base::getIdentityOnTxM(RefMat out, const ConstRefVec& x) const
  {
    mnf_assert(out.rows() == tangentDim_);
    mnf_assert(out.cols() == tangentDim_);
    mnf_assert(x.size() == representationDim());
    getIdentityOnTxM_(out, x);
  }

  void Manifold_Base::limitMap(RefVec out) const
  {
    mnf_assert(out.size() == tangentDim_);
    limitMap_(out);
  }

  void Manifold_Base::getTypicalMagnitude(RefVec out) const
  {
    mnf_assert(out.size() == tangentDim_);
    out = typicalMagnitude_;
  }

  Eigen::VectorXd Manifold_Base::getTypicalMagnitude() const
  {
    Eigen::VectorXd out(tangentDim_);
    getTypicalMagnitude(out);
    return out;
  }

  void Manifold_Base::setTypicalMagnitude(const double& magnitude)
  {
    Eigen::VectorXd temp = Eigen::VectorXd::Constant(tangentDim(), magnitude);
    setTypicalMagnitude(temp);
  }

  void Manifold_Base::setTypicalMagnitude(const ConstRefVec& out)
  {
    typicalMagnitude_ = out;
  }

  void Manifold_Base::getTrustMagnitude(RefVec out) const
  {
    mnf_assert(out.size() == tangentDim_);
    out = trustMagnitude_;
  }

  Eigen::VectorXd Manifold_Base::getTrustMagnitude() const
  {
    Eigen::VectorXd out(tangentDim_);
    getTrustMagnitude(out);
    return out;
  }

  void Manifold_Base::setTrustMagnitude(const double& magnitude)
  {
    Eigen::VectorXd temp = Eigen::VectorXd::Constant(tangentDim(), magnitude);
    setTrustMagnitude(temp);
  }

  void Manifold_Base::setTrustMagnitude(const ConstRefVec& out)
  {
    trustMagnitude_ = out;
  }

  void Manifold_Base::lock() const
  {
    lock_ = true;
  }

  void Manifold_Base::testLock() const
  {
    if (lock_)
      throw std::runtime_error("Either a point or a compound manifold is relying on this manifold, you can't modify it anymore.");
  }

  bool Manifold_Base::isLocked() const
  {
    return lock_;
  }

  const std::string& Manifold_Base::name() const { return name_;}

  std::string& Manifold_Base::name() { return name_;}

  long Manifold_Base::getInstanceId() const
  {
    return this->instanceId_;
  }

  bool Manifold_Base::isSameType(const Manifold_Base& other) const
  {
    return this->getTypeId() == other.getTypeId();
  }

}

