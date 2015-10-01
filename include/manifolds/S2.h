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

#pragma once
#include <stdexcept>

#include <manifolds/defs.h>
#include <manifolds/Manifold.h>
#include <manifolds/utils.h>

namespace mnf
{
  class S2_Base;

  /// \brief Manifold representing the 3-dimensional Sphere, also
  /// known as S2.
  /// All the equations in this class are provided by Manopt
  class MANIFOLDS_API S2 : public Manifold
  {
  public:
    S2(std::shared_ptr<S2_Base> p);
    S2();
    S2(double mag);
    S2(const ConstRefVec& mag);
    virtual S2 copy() const;

    void logarithm (RefVec out, const ConstRefVec& x, const ConstRefVec& y) const;
    double distance (const ConstRefVec& x, const ConstRefVec& y) const;
    /// \brief projects each row of \ in on TxM
    void projRows (RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    /// \brief projects each cols of \ in on TxM
    void projCols (RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    /// \brief projects vector \ in on TxM
    void projVec (RefVec out, const ConstRefVec& in, const ConstRefVec& x) const;
    void rand(RefVec out) const;
    void randVec(RefVec out, const ConstRefVec& x) const;
    Eigen::Vector3d randVec(const ConstRefVec& x) const;

  };
}
