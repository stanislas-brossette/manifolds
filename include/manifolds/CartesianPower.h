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

#include <manifolds/defs.h>
#include <manifolds/CartesianProduct.h>

namespace mnf
{
  /// \brief Manifold representing the cartesian product of n times the same manifold
  class MANIFOLDS_API CartesianPower : public CartesianProduct
  {
  public:
    /// \brief Constructor of the \f$ M^n m2\f$
    CartesianPower(const Manifold& M, const int n);
  };

}


