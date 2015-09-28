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
    : Manifold(2, 3, 3)
  {
  }

  S2::S2(double magnitude)
    : Manifold(2, 3, 3)
  {
  }

  S2::S2(const ConstRefVec& magnitude)
    : Manifold(2, 3, 3)
  {
  }

  size_t S2::numberOfSubmanifolds() const
  {
  }

  bool S2::isElementary() const
  {
  }

  const Manifold& S2::operator()(size_t i) const
  {
  }

  std::string S2::toString(const ConstRefVec& val, const std::string& prefix, const Eigen::IOFormat& fmt) const
  {
  }

  void S2::createRandomPoint_(RefVec out, double ) const
  {
  }

  void S2::logarithm (RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
  }

  double S2::distance (const ConstRefVec& x, const ConstRefVec& y) const
  {
  }

  void S2::projRows(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
  }
  void S2::projCols(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
  }
  void S2::projVec(RefVec out, const ConstRefVec& in, const ConstRefVec& x) const
  {
  }

  void S2::rand(RefVec out) const
  {
  }

  void S2::randVec(RefVec out, const ConstRefVec& x) const
  {
  }
  Eigen::Vector3d S2::randVec(const ConstRefVec& x) const
  {
  }

  void S2::getTypicalMagnitude_(RefVec out) const
  {
  }

  void S2::setTypicalMagnitude(double magnitude)
  {
  }

  void S2::setTypicalMagnitude(const ConstRefVec& out)
  {
  }

  void S2::getTrustMagnitude_(RefVec out) const
  {
  }

  void S2::setTrustMagnitude(const double& magnitude)
  {
  }

  void S2::setTrustMagnitude(const ConstRefVec& out)
  {
  }

  long S2::getTypeId() const
  {
  }
}
