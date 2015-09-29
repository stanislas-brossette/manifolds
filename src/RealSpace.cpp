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
    : Manifold(std::shared_ptr<Manifold_Base>(new RealSpace_Base(n)))
  {
  }

  RealSpace::RealSpace(Index n, double magnitude)
    : Manifold(std::shared_ptr<Manifold_Base>(new RealSpace_Base(n, magnitude)))
  {
  }
  RealSpace::RealSpace(Index n, const ConstRefVec& magnitude)
    : Manifold(std::shared_ptr<Manifold_Base>(new RealSpace_Base(n, magnitude)))
  {
  }

  size_t RealSpace::numberOfSubmanifolds() const
  {
    return std::static_pointer_cast<RealSpace_Base>(manifoldBase_)->numberOfSubmanifolds();
  }

  bool RealSpace::isElementary() const
  {
    return std::static_pointer_cast<RealSpace_Base>(manifoldBase_)->isElementary();
  }

  Manifold RealSpace::operator()(size_t i) const
  {
    return makeManifold(std::const_pointer_cast<Manifold_Base>(manifoldBase_->operator()(i)));;
  }

  std::string RealSpace::toString(const ConstRefVec& val, const std::string& prefix, const Eigen::IOFormat& fmt) const
  {
    return std::static_pointer_cast<RealSpace_Base>(manifoldBase_)->toString(val, prefix, fmt);
  }

  void RealSpace::getTypicalMagnitude_(RefVec out) const
  {
    return std::static_pointer_cast<RealSpace_Base>(manifoldBase_)->getTypicalMagnitude_(out);
  }

  void RealSpace::setTypicalMagnitude(double magnitude)
  {
    return std::static_pointer_cast<RealSpace_Base>(manifoldBase_)->setTypicalMagnitude(magnitude);
  }

  void RealSpace::setTypicalMagnitude(const ConstRefVec& out)
  {
    return std::static_pointer_cast<RealSpace_Base>(manifoldBase_)->setTypicalMagnitude(out);
  }

  void RealSpace::getTrustMagnitude_(RefVec out) const
  {
    return std::static_pointer_cast<RealSpace_Base>(manifoldBase_)->getTrustMagnitude_(out);
  }

  void RealSpace::setTrustMagnitude(const double& magnitude)
  {
    return std::static_pointer_cast<RealSpace_Base>(manifoldBase_)->setTrustMagnitude(magnitude);
  }

  void RealSpace::setTrustMagnitude(const ConstRefVec& out)
  {
    return std::static_pointer_cast<RealSpace_Base>(manifoldBase_)->setTrustMagnitude(out);
  }

  long RealSpace::getTypeId() const
  {
    return std::static_pointer_cast<RealSpace_Base>(manifoldBase_)->getTypeId();
  }
}
