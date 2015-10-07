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
    : Manifold(std::shared_ptr<internal::Manifold_Base>(new internal::CartesianProduct_Base()))
{std::cout << "CartesianProduct::CartesianProduct()" << std::endl;
}

CartesianProduct::CartesianProduct(CartesianProduct&& m)
    : Manifold(std::move(m))
{std::cout << "CartesianProduct::CartesianProduct(CartesianProduct&& m)" << std::endl;
}

CartesianProduct::CartesianProduct(Manifold&& m)
    : Manifold(nullptr)
{std::cout << "CartesianProduct::CartesianProduct(Manifold&& m)" << std::endl;
  CartesianProduct(dynamic_cast<CartesianProduct&&>(m));
  if (dynamic_cast<internal::CartesianProduct_Base*>(m.getNonConstPtr().get()) != nullptr)
  {
    ptr_ = m.getNonConstPtr();
    m.getNonConstPtr() = nullptr;
  }
  else
    throw std::runtime_error("Up");
}

CartesianProduct::CartesianProduct(std::initializer_list<Manifold*> listM)
    : Manifold(std::shared_ptr<internal::Manifold_Base>(new internal::CartesianProduct_Base()))
{std::cout << "CartesianProduct::CartesianProduct(std::initializer_list<Manifold*> listM)" << std::endl;
  for (auto m : listM) 
    multiply(*m);
}

CartesianProduct::CartesianProduct(Manifold& m1, Manifold& m2)
    : Manifold(std::shared_ptr<internal::Manifold_Base>(new internal::CartesianProduct_Base()))
{std::cout << "CartesianProduct::CartesianProduct(Manifold& m1, Manifold& m2)" << std::endl;
  multiply(m1);
  multiply(m2);
}

Manifold CartesianProduct::shallowCopy()
{std::cout << "Manifold CartesianProduct::shallowCopy()" << std::endl;
  CartesianProduct m;
  for(std::size_t i=0; i < subManifolds_.size(); ++i)
    m.multiply(subManifolds_[i].shallowCopy());
  return std::move(m);
}

Manifold CartesianProduct::deepCopy() const
{std::cout << "Manifold CartesianProduct::deepCopy()" << std::endl;
  CartesianProduct m;
  for(std::size_t i=0; i < subManifolds_.size(); ++i)
    m.multiply(subManifolds_[i].deepCopy());
  return std::move(m);
}

CartesianProduct& CartesianProduct::multiply(Manifold& m)
{std::cout << "CartesianProduct& CartesianProduct::multiply(Manifold& m)" << std::endl;
  std::static_pointer_cast<internal::CartesianProduct_Base>(ptr_)->multiply(m.getNonConstPtr());
  subManifolds_.push_back(m.shallowCopy());
  return *this;
}

CartesianProduct& CartesianProduct::multiply(Manifold&& m)
{std::cout << "CartesianProduct& CartesianProduct::multiply(Manifold&& m)" << std::endl;
  std::static_pointer_cast<internal::CartesianProduct_Base>(ptr_)
      ->multiply(m.getNonConstPtr());
  subManifolds_.push_back(std::move(m));
  return *this;
}

const Manifold& CartesianProduct::operator()(size_t i) const
{std::cout << "const Manifold& CartesianProduct::operator()" << std::endl;
  return subManifolds_[i];
}

Manifold& CartesianProduct::operator()(size_t i)
{std::cout << "Manifold& CartesianProduct::operator()" << std::endl;
  return subManifolds_[i];
}

//CartesianProduct CartesianProduct::copy() const
//{
  //return CartesianProduct(std::make_shared<internal::CartesianProduct_Base>(
      //std::static_pointer_cast<internal::CartesianProduct_Base>(ptr_)->copy()));
//}
}
