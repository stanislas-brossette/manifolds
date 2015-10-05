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

#define _USE_MATH_DEFINES
#include <math.h>

#include <manifolds/defs.h>
#include <manifolds/Manifold.h>
#include <manifolds/ReusableTemporaryMap.h>
#include <manifolds/utils.h>
#include <manifolds/SO3_Base.h>

namespace mnf
{
/// \brief Manifold representing the space of 3-dimensional rotations, also
/// known as SO(3). It is templated by its map
template <typename Map>
class SO3 : public Manifold
{
 public:
  SO3(std::shared_ptr<SO3_Base<Map>> p);
  SO3();
  SO3(double magnitude);
  SO3(const ConstRefVec& magnitude);
  virtual SO3 copy() const;
};

// Implementations of the methods
template <typename Map>
inline SO3<Map>::SO3(std::shared_ptr<SO3_Base<Map>> p)
    : Manifold(std::static_pointer_cast<SO3_Base<Map>>(p))
{
}
template <typename Map>
inline SO3<Map>::SO3()
    : Manifold(std::shared_ptr<Manifold_Base>(new SO3_Base<Map>()))
{
}
template <typename Map>
inline SO3<Map>::SO3(double magnitude)
    : Manifold(std::shared_ptr<Manifold_Base>(new SO3_Base<Map>(magnitude)))
{
}
template <typename Map>
inline SO3<Map>::SO3(const ConstRefVec& magnitude)
    : Manifold(std::shared_ptr<Manifold_Base>(new SO3_Base<Map>(magnitude)))
{
}
template <typename Map>
inline SO3<Map> SO3<Map>::copy() const
{
  return SO3<Map>(std::make_shared<SO3_Base<Map>>(
      std::static_pointer_cast<SO3_Base<Map>>(ptr_)->copy()));
}
}

