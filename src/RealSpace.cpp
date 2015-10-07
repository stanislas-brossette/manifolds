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

#include <manifolds/RealSpace.h>
#include <manifolds/mnf_assert.h>

namespace mnf
{
RealSpace::RealSpace(std::shared_ptr<internal::RealSpace_Base> p)
    : Manifold(std::static_pointer_cast<internal::Manifold_Base>(p))
{
  std::cout << "RealSpace::RealSpace(std::shared_ptr<internal::RealSpace_Base> p)" << std::endl;
}

RealSpace::RealSpace(RealSpace&& m)
    : Manifold(std::move(m))
{
  std::cout << "RealSpace::RealSpace(RealSpace&& m)" << std::endl;
}

RealSpace::RealSpace(Manifold&& m)
    : Manifold(nullptr)
{
  std::cout << "RealSpace::RealSpace(Manifold&& m)" << std::endl;
  if (dynamic_cast<internal::RealSpace_Base*>(m.getNonConstPtr().get()) != nullptr)
  {
    ptr_ = m.getNonConstPtr();
    m.getNonConstPtr() = nullptr;
  }
  else
    throw std::runtime_error("Up");
}

RealSpace::RealSpace(Index n)
    : Manifold(std::shared_ptr<internal::Manifold_Base>(new internal::RealSpace_Base(n)))
{
  std::cout << "RealSpace::RealSpace(Index n)" << std::endl;
}

RealSpace::RealSpace(Index n, double magnitude)
    : Manifold(std::shared_ptr<internal::Manifold_Base>(new internal::RealSpace_Base(n, magnitude)))
{
  std::cout << "RealSpace::RealSpace(Index n, double mag)" << std::endl;
}
RealSpace::RealSpace(Index n, const ConstRefVec& magnitude)
    : Manifold(std::shared_ptr<internal::Manifold_Base>(new internal::RealSpace_Base(n, magnitude)))
{
  std::cout << "RealSpace::RealSpace(Index n, Vec mag)" << std::endl;
}

//Manifold RealSpace::shallowCopy()
//{
  //std::cout << "RealSpace::shallowCopy()" << std::endl;
  //return (RealSpace(ptr_));
//}

//Manifold RealSpace::deepCopy() const
//{
  //std::cout << "RealSpace::deepCopy()" << std::endl;
  //return (RealSpace(ptr_->clone()));
//}

RealSpace RealSpace::copy() const
{
  std::cout << "RealSpace RealSpace::copy() const" << std::endl;
  return RealSpace(std::make_shared<internal::RealSpace_Base>(
      std::static_pointer_cast<internal::RealSpace_Base>(ptr_)->copy()));
}
}
