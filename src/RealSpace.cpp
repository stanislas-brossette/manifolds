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
  }

  RealSpace::RealSpace(Index n, double magnitude)
    : Manifold(n, n, n)
  {
  }
  RealSpace::RealSpace(Index n, const ConstRefVec& magnitude)
    : Manifold(n, n, n)
  {
  }

  bool RealSpace::isInM_(const Eigen::VectorXd& val , double ) const
  {
  }

  void RealSpace::forceOnM_(RefVec out, const ConstRefVec& in) const
  {
  }

  void RealSpace::getIdentityOnTxM_(RefMat out, const ConstRefVec&) const
  {
  }

  size_t RealSpace::numberOfSubmanifolds() const
  {
  }

  bool RealSpace::isElementary() const
  {
  }

  const Manifold& RealSpace::operator()(size_t i) const
  {
  }

  std::string RealSpace::toString(const ConstRefVec& val, const std::string& prefix, const Eigen::IOFormat& fmt) const
  {
  }

  void RealSpace::createRandomPoint_(RefVec out, double coeff) const
  {
  }

  void RealSpace::retractation_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
  }

  void RealSpace::pseudoLog_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
  }

  void RealSpace::pseudoLog0_(RefVec out, const ConstRefVec& x) const
  {
  }

  void RealSpace::setZero_(RefVec out) const
  {
  }

  Eigen::MatrixXd RealSpace::diffRetractation_(const ConstRefVec& ) const
  {
  }

  void RealSpace::applyDiffRetractation_(RefMat out, const ConstRefMat& in, const ConstRefVec& ) const
  {
  }

  Eigen::MatrixXd RealSpace::diffPseudoLog0_(const ConstRefVec&) const
  {
  }

  void RealSpace::applyDiffPseudoLog0_(RefMat out, const ConstRefMat& in, const ConstRefVec& ) const
  {
  }

  void RealSpace::applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec& ) const
  {
  }

  void RealSpace::applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec& ) const
  {
  }

  void RealSpace::applyInvTransportOnTheRight_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec&) const
  {
  }

  void RealSpace::tangentConstraint_(RefMat, const ConstRefVec&) const
  {
  }

  bool RealSpace::isInTxM_(const ConstRefVec&, const ConstRefVec&, const double& ) const
  {
  }

  void RealSpace::forceOnTxM_(RefVec out, const ConstRefVec& in, const ConstRefVec&) const
  {
  }

  void RealSpace::limitMap_(RefVec out) const
  {
  }

  void RealSpace::getTypicalMagnitude_(RefVec out) const
  {
  }

  void RealSpace::setTypicalMagnitude(double magnitude)
  {
  }

  void RealSpace::setTypicalMagnitude(const ConstRefVec& out)
  {
  }

  void RealSpace::getTrustMagnitude_(RefVec out) const
  {
  }

  void RealSpace::setTrustMagnitude(const double& magnitude)
  {
  }

  void RealSpace::setTrustMagnitude(const ConstRefVec& out)
  {
  }

  long RealSpace::getTypeId() const
  {
  }

  Manifold_ptr RealSpace::getNewCopy() const
  {
  }
}
