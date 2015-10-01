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

#include <limits>

#include <manifolds/RealSpace_Base.h>
#include <manifolds/mnf_assert.h>

namespace mnf
{
  RealSpace_Base::RealSpace_Base(Index n)
    : Manifold_Base(n, n, n)
  {
    name() = "R" + std::to_string( n );
    typicalMagnitude_.resize(n);
    trustMagnitude_.resize(n);
    Eigen::VectorXd temp = Eigen::VectorXd::Constant(n, 1.0);
    setTypicalMagnitude(temp);
    setTrustMagnitude(temp);
  }

  RealSpace_Base::RealSpace_Base(Index n, double magnitude)
    : Manifold_Base(n, n, n)
  {
    name() = "R" + std::to_string( n );
    typicalMagnitude_.resize(n);
    trustMagnitude_.resize(n);
    Eigen::VectorXd temp = Eigen::VectorXd::Constant(n, magnitude);
    setTypicalMagnitude(temp);
    setTrustMagnitude(temp);
  }
  RealSpace_Base::RealSpace_Base(Index n, const ConstRefVec& magnitude)
    : Manifold_Base(n, n, n)
  {
    mnf_assert(magnitude.size() == n && "magnitude on R^n must be of size n");
    name() = "R" + std::to_string( n );
    typicalMagnitude_.resize(n);
    setTypicalMagnitude (magnitude);
    trustMagnitude_.resize(n);
    setTrustMagnitude (magnitude);
  }

  RealSpace_Base RealSpace_Base::copy() const
  {
    return RealSpace_Base(*this);
  }

  bool RealSpace_Base::isInM_(const Eigen::VectorXd& val , double ) const
  {
    bool out( dim() == val.size());
    return out;
  }

  void RealSpace_Base::forceOnM_(RefVec out, const ConstRefVec& in) const
  {
    out = in;
  }

  void RealSpace_Base::getIdentityOnTxM_(RefMat out, const ConstRefVec&) const
  {
    out.setIdentity();
  }

  size_t RealSpace_Base::numberOfSubmanifolds() const
  {
    return 1;
  }

  bool RealSpace_Base::isElementary() const
  {
    return true;
  }

  std::shared_ptr<const Manifold_Base> RealSpace_Base::operator()(size_t ) const
  {
    return shared_from_this();
  }

  std::string RealSpace_Base::toString(const ConstRefVec& val, const std::string& prefix, const Eigen::IOFormat& fmt) const
  {
    std::stringstream ss;
    ss << prefix << val.transpose().format(fmt);
    return ss.str();
  }

  void RealSpace_Base::createRandomPoint_(RefVec out, double coeff) const
  {
    out = coeff*Eigen::VectorXd::Random(representationDim());
  }

  void RealSpace_Base::retractation_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    out = x + v;
  }

  void RealSpace_Base::pseudoLog_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    out = y - x;
  }

  void RealSpace_Base::pseudoLog0_(RefVec out, const ConstRefVec& x) const
  {
    out = x;
  }

  void RealSpace_Base::setZero_(RefVec out) const
  {
    out.setZero();
  }

  Eigen::MatrixXd RealSpace_Base::diffRetractation_(const ConstRefVec& ) const
  {
    return Eigen::MatrixXd::Identity(representationDim(),dim());
  }

  void RealSpace_Base::applyDiffRetractation_(RefMat out, const ConstRefMat& in, const ConstRefVec& ) const
  {
    out = in;
  }

  Eigen::MatrixXd RealSpace_Base::diffPseudoLog0_(const ConstRefVec&) const
  {
    return Eigen::MatrixXd::Identity(representationDim(),dim());
  }

  void RealSpace_Base::applyDiffPseudoLog0_(RefMat out, const ConstRefMat& in, const ConstRefVec& ) const
  {
    out = in;
  }

  void RealSpace_Base::applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec& ) const
  {
    out = in;
  }

  void RealSpace_Base::applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec& ) const
  {
    out = in;
  }

  void RealSpace_Base::applyInvTransportOnTheRight_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec&) const
  {
    out = in;
  }

  void RealSpace_Base::tangentConstraint_(RefMat, const ConstRefVec&) const
  {
    //matrix is 0xt, no need to fill it.
  }

  bool RealSpace_Base::isInTxM_(const ConstRefVec&, const ConstRefVec&, const double& ) const
  {
    return true;
  }

  void RealSpace_Base::forceOnTxM_(RefVec out, const ConstRefVec& in, const ConstRefVec&) const
  {
    out = in;
  }

  void RealSpace_Base::limitMap_(RefVec out) const
  {
    out.setConstant(std::numeric_limits<double>::infinity());
  }

  long RealSpace_Base::getTypeId() const
  {
    long typeId = utils::hash::computeHash("RealSpace");
    return typeId;
  }

  Manifold_Base_ptr RealSpace_Base::getNewCopy() const
  {
    RealSpace_Base_ptr copy(new RealSpace_Base(*this));
    copy->instanceId_ = this->instanceId_;

    return copy;
  }
}
