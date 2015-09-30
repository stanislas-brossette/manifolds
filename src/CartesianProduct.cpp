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

  //CartesianProduct::CartesianProduct(const std::initializer_list<Manifold_Base*> m)
    //: Manifold(std::shared_ptr<Manifold_Base>(new CartesianProduct_Base(m)))
  //{
    //std::cout << "\n\n\n Need to fix:\n CartesianProduct::CartesianProduct(const std::initializer_list<Manifold_Base*> m) should be CartesianProduct::CartesianProduct(const std::initializer_list<Manifold*> m)\n\n\n" << std::endl;
  //}

  CartesianProduct::CartesianProduct(const Manifold& m1, const Manifold& m2)
    : Manifold(std::shared_ptr<Manifold_Base>(new CartesianProduct_Base(m1, m2)))
  {
  }

  //void CartesianProduct::getIdentityOnTxM_(RefMat out, const ConstRefVec& x) const
  //{
    //return std::static_pointer_cast<CartesianProduct_Base>(manifoldBase_)->getIdentityOnTxM_(out, x);
  //}

  CartesianProduct& CartesianProduct::multiply(const Manifold& m)
  {
    std::static_pointer_cast<CartesianProduct_Base>(manifoldBase_)->multiply(m);
    return *this;
    //std::static_pointer_cast<CartesianProduct_Base>(manifoldBase_)->multiply(m.getManifoldBase());
  }

  //size_t CartesianProduct::numberOfSubmanifolds() const
  //{
    //return std::static_pointer_cast<CartesianProduct_Base>(manifoldBase_)->numberOfSubmanifolds();
  //}

  //bool CartesianProduct::isElementary() const
  //{
    //return std::static_pointer_cast<CartesianProduct_Base>(manifoldBase_)->isElementary();
  //}

  //void CartesianProduct::display(const std::string& prefix) const
  //{
    //return std::static_pointer_cast<CartesianProduct_Base>(manifoldBase_)->display(prefix);
  //}

  //const Manifold& CartesianProduct::operator()(const size_t i) const
  //{
    //return makeManifold(std::const_pointer_cast<CartesianProduct_Base>(manifoldBase_->operator()(i)));;
  //}

  //std::string CartesianProduct::toString(const ConstRefVec& val, const std::string& prefix, const Eigen::IOFormat& fmt) const
  //{
    //return std::static_pointer_cast<CartesianProduct_Base>(manifoldBase_)->toString(val, prefix, fmt);
  //}
}
