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

#include <manifolds/defs.h>
#include <manifolds/CartesianProduct.h>
#include <manifolds/CartesianProduct_Base.h>
#include <manifolds/mnf_assert.h>

namespace mnf
{
  CartesianProduct::CartesianProduct()
    : Manifold(std::shared_ptr<Manifold_Base>(new CartesianProduct_Base()))
  {
  }

  CartesianProduct::CartesianProduct(const std::initializer_list<Manifold*> m)
    : Manifold(std::shared_ptr<Manifold_Base>(new CartesianProduct_Base()))
  {
    for(auto s:m)
      multiply(s->getNonConstPtr());
  }

  CartesianProduct::CartesianProduct(const Manifold& m1, const Manifold& m2)
    : Manifold(std::shared_ptr<Manifold_Base>(
          new CartesianProduct_Base(m1.getNonConstPtr(), m2.getNonConstPtr())))
  {
  }

  CartesianProduct& CartesianProduct::multiply(const Manifold& m)
  {
    std::static_pointer_cast<CartesianProduct_Base>(ptr_)->multiply(m.getNonConstPtr());
    return *this;
  }
}
