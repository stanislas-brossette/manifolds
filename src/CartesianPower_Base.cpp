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

#include <manifolds/CartesianPower_Base.h>

namespace mnf
{
  CartesianPower_Base::CartesianPower_Base(std::shared_ptr<Manifold_Base> M, const int n)
    : CartesianProduct_Base()
  {
    for (int i = 0; i < n; ++i)
      this->multiply(std::make_shared<const Manifold_Base>(*M.get()));
  }
}


