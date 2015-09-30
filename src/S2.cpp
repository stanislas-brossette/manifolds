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
#define _USE_MATH_DEFINES
#include <math.h>

#include <manifolds/utils.h>
#include <manifolds/S2.h>
#include <manifolds/S2_Base.h>
#include <manifolds/mnf_assert.h>

namespace mnf
{
  S2::S2()
    : Manifold(std::shared_ptr<Manifold_Base>(new S2_Base()))
  {
  }

  S2::S2(double magnitude)
    : Manifold(std::shared_ptr<Manifold_Base>(new S2_Base(magnitude)))
  {
  }

  S2::S2(const ConstRefVec& magnitude)
    : Manifold(std::shared_ptr<Manifold_Base>(new S2_Base(magnitude)))
  {
  }

  size_t S2::numberOfSubmanifolds() const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->numberOfSubmanifolds();
  }

  bool S2::isElementary() const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->isElementary();
  }

  Manifold S2::operator()(size_t i) const
  {
    return makeManifold(std::const_pointer_cast<Manifold_Base>(manifoldBase_->operator()(i)));;
  }

  std::string S2::toString(const ConstRefVec& val, const std::string& prefix, const Eigen::IOFormat& fmt) const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->toString(val, prefix, fmt);
  }

  void S2::createRandomPoint_(RefVec out, double d) const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->createRandomPoint_(out, d);
  }

  void S2::logarithm (RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->logarithm (out, x, y);
  }

  double S2::distance (const ConstRefVec& x, const ConstRefVec& y) const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->distance (x, y);
  }

  void S2::projRows(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->projRows(out, in, x);
  }
  void S2::projCols(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->projCols(out, in, x);
  }
  void S2::projVec(RefVec out, const ConstRefVec& in, const ConstRefVec& x) const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->projVec(out, in, x);
  }

  void S2::rand(RefVec out) const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->rand(out);
  }

  void S2::randVec(RefVec out, const ConstRefVec& x) const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->randVec(out, x);
  }
  Eigen::Vector3d S2::randVec(const ConstRefVec& x) const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->randVec(x);
  }

  void S2::getTypicalMagnitude_(RefVec out) const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->getTypicalMagnitude_(out);
  }

  void S2::setTypicalMagnitude(double magnitude)
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->setTypicalMagnitude(magnitude);
  }

  void S2::setTypicalMagnitude(const ConstRefVec& out)
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->setTypicalMagnitude(out);
  }

  void S2::getTrustMagnitude_(RefVec out) const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->getTrustMagnitude_(out);
  }

  void S2::setTrustMagnitude(const double& magnitude)
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->setTrustMagnitude(magnitude);
  }

  void S2::setTrustMagnitude(const ConstRefVec& out)
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->setTrustMagnitude(out);
  }

  long S2::getTypeId() const
  {
    return std::static_pointer_cast<S2_Base>(manifoldBase_)->getTypeId();
  }
}
