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
#include <manifolds/mnf_assert.h>

namespace mnf
{
  CartesianProduct::CartesianProduct()
    : Manifold(std::make_shared<CartesianProduct_Base>())
  {
  }

  CartesianProduct::CartesianProduct(const std::initializer_list<Manifold*> m)
    : Manifold(std::make_shared<CartesianProduct_Base>(m))
  {
  }

  CartesianProduct::CartesianProduct(const Manifold& m1, const Manifold& m2)
    : Manifold(std::make_shared<CartesianProduct_Base>(m1, m2))
  {
  }

  void CartesianProduct::getIdentityOnTxM_(RefMat out, const ConstRefVec& x) const
  {
    return std::static_pointer_cast<CartesianProduct_Base>(manifoldBase_)->getIdentityOnTxM_(out, x);
  }

  CartesianProduct& CartesianProduct::multiply(const Manifold& m)
  {
    return std::static_pointer_cast<CartesianProduct_Base>(manifoldBase_)->multiply(m);
  }

  size_t CartesianProduct::numberOfSubmanifolds() const
  {
    return std::static_pointer_cast<CartesianProduct_Base>(manifoldBase_)->numberOfSubmanifolds();
  }

  bool CartesianProduct::isElementary() const
  {
    return std::static_pointer_cast<CartesianProduct_Base>(manifoldBase_)->isElementary();
  }

  void CartesianProduct::display(const std::string& prefix) const
  {
    return std::static_pointer_cast<CartesianProduct_Base>(manifoldBase_)->display(prefix);
  }

  const Manifold& CartesianProduct::operator()(const size_t i) const
  {
    return std::static_pointer_cast<CartesianProduct_Base>(manifoldBase_)->operator();
  }

  std::string CartesianProduct::toString(const ConstRefVec& val, const std::string& prefix, const Eigen::IOFormat& fmt) const
  {
    return std::static_pointer_cast<CartesianProduct_Base>(manifoldBase_)->toString(val, prefix, fmt);
  }
}
