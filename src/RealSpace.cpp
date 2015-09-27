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

#include <manifolds/RealSpace.h>
#include <manifolds/mnf_assert.h>

namespace mnf
{
  RealSpace::RealSpace(Index n)
    : Manifold(n, n, n)
  {
    name() = "R" + std::to_string( n );
    typicalMagnitude_.resize(n);
    trustMagnitude_.resize(n);
    Eigen::VectorXd temp = Eigen::VectorXd::Constant(n, 1.0);
    setTypicalMagnitude(temp);
    setTrustMagnitude(temp);
  }

  RealSpace::RealSpace(Index n, double magnitude)
    : Manifold(n, n, n)
  {
    name() = "R" + std::to_string( n );
    typicalMagnitude_.resize(n);
    trustMagnitude_.resize(n);
    Eigen::VectorXd temp = Eigen::VectorXd::Constant(n, magnitude);
    setTypicalMagnitude(temp);
    setTrustMagnitude(temp);
  }
  RealSpace::RealSpace(Index n, const ConstRefVec& magnitude)
    : Manifold(n, n, n)
  {
    mnf_assert(magnitude.size() == n && "magnitude on R^n must be of size n");
    name() = "R" + std::to_string( n );
    typicalMagnitude_.resize(n);
    setTypicalMagnitude (magnitude);
    trustMagnitude_.resize(n);
    setTrustMagnitude (magnitude);
  }

  bool RealSpace::isInM_(const Eigen::VectorXd& val , double ) const
  {
    bool out( dim() == val.size());
    return out;
  }

  void RealSpace::forceOnM_(RefVec out, const ConstRefVec& in) const
  {
    out = in;
  }

  void RealSpace::getIdentityOnTxM_(RefMat out, const ConstRefVec&) const
  {
    out.setIdentity();
  }

  size_t RealSpace::numberOfSubmanifolds() const
  {
    return 1;
  }

  bool RealSpace::isElementary() const
  {
    return true;
  }

  const Manifold& RealSpace::operator()(size_t i) const
  {
    mnf_assert(i < 1 && "invalid index");
    return *this;
  }

  std::string RealSpace::toString(const ConstRefVec& val, const std::string& prefix, const Eigen::IOFormat& fmt) const
  {
    std::stringstream ss;
    ss << prefix << val.transpose().format(fmt);
    return ss.str();
  }

  void RealSpace::createRandomPoint_(RefVec out, double coeff) const
  {
    out = coeff*Eigen::VectorXd::Random(representationDim());
  }

  void RealSpace::retractation_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    out = x + v;
  }

  void RealSpace::pseudoLog_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    out = y - x;
  }

  void RealSpace::pseudoLog0_(RefVec out, const ConstRefVec& x) const
  {
    out = x;
  }

  void RealSpace::setZero_(RefVec out) const
  {
    out.setZero();
  }

  Eigen::MatrixXd RealSpace::diffRetractation_(const ConstRefVec& ) const
  {
    return Eigen::MatrixXd::Identity(representationDim(),dim());
  }

  void RealSpace::applyDiffRetractation_(RefMat out, const ConstRefMat& in, const ConstRefVec& ) const
  {
    out = in;
  }

  Eigen::MatrixXd RealSpace::diffPseudoLog0_(const ConstRefVec&) const
  {
    return Eigen::MatrixXd::Identity(representationDim(),dim());
  }

  void RealSpace::applyDiffPseudoLog0_(RefMat out, const ConstRefMat& in, const ConstRefVec& ) const
  {
    out = in;
  }

  void RealSpace::applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec& ) const
  {
    out = in;
  }

  void RealSpace::applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec& ) const
  {
    out = in;
  }

  void RealSpace::applyInvTransportOnTheRight_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec&) const
  {
    out = in;
  }

  void RealSpace::tangentConstraint_(RefMat, const ConstRefVec&) const
  {
    //matrix is 0xt, no need to fill it.
  }

  bool RealSpace::isInTxM_(const ConstRefVec&, const ConstRefVec&, const double& ) const
  {
    return true;
  }

  void RealSpace::forceOnTxM_(RefVec out, const ConstRefVec& in, const ConstRefVec&) const
  {
    out = in;
  }

  void RealSpace::limitMap_(RefVec out) const
  {
    out.setConstant(std::numeric_limits<double>::infinity());
  }

  void RealSpace::getTypicalMagnitude_(RefVec out) const
  {
    out = typicalMagnitude_;
  }

  void RealSpace::setTypicalMagnitude(double magnitude)
  {
    Eigen::VectorXd temp = Eigen::VectorXd::Constant(tangentDim(), magnitude);
    setTypicalMagnitude(temp);
  }

  void RealSpace::setTypicalMagnitude(const ConstRefVec& out)
  {
    testLock();
    typicalMagnitude_ = out;
  }

  void RealSpace::getTrustMagnitude_(RefVec out) const
  {
    out = trustMagnitude_;
  }

  void RealSpace::setTrustMagnitude(const double& magnitude)
  {
    Eigen::VectorXd temp = Eigen::VectorXd::Constant(tangentDim(), magnitude);
    setTrustMagnitude(temp);
  }

  void RealSpace::setTrustMagnitude(const ConstRefVec& out)
  {
    testLock();
    trustMagnitude_ = out;
  }

  long RealSpace::getTypeId() const
  {
    long typeId = utils::hash::computeHash("RealSpace");
    return typeId;
  }

  Manifold* RealSpace::getNewCopy() const
  {
    RealSpace* copy = new RealSpace(*this);
    copy->instanceId_ = this->instanceId_;

    return copy;
  }
}
