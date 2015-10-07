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

#include <iostream>
#include <manifolds/defs.h>
#include <manifolds/Manifold.h>
#include <manifolds/RealSpace_Base.h>
#include <manifolds/utils.h>

namespace mnf
{
/// \brief Manifold representing the space of real numbers of dimension n
/// \f$\mathbb{R}^n\f$
class MANIFOLDS_API RealSpace : public Manifold
{
 private:
  RealSpace(std::shared_ptr<internal::RealSpace_Base> p);

 public:
  RealSpace(RealSpace&& m);
  RealSpace(Manifold&& m);

  //virtual Manifold shallowCopy();
  //virtual Manifold deepCopy() const;

  /// \brief Constructor
  /// \param n the dimension of the realspace \f$\mathbb{R}^n\f$
  RealSpace(Index n);
  RealSpace(Index n, double magnitude);
  RealSpace(Index n, const ConstRefVec& magnitude);

  /// @brief Returns a deep copy of this manifold.
  /// The output contains a copy of this Manifold_Base that lives independently
  /// from this manifold
  RealSpace copy() const;
};
}

