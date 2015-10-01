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
  RealSpace::RealSpace(std::shared_ptr<RealSpace_Base> p)
    : Manifold(std::static_pointer_cast<Manifold_Base>(p))
  {
  }

  RealSpace::RealSpace(Index n)
    : Manifold(std::shared_ptr<Manifold_Base>(new RealSpace_Base(n)))
  {
  }

  RealSpace::RealSpace(Index n, double magnitude)
    : Manifold(std::shared_ptr<Manifold_Base>(new RealSpace_Base(n, magnitude)))
  {
  }
  RealSpace::RealSpace(Index n, const ConstRefVec& magnitude)
    : Manifold(std::shared_ptr<Manifold_Base>(new RealSpace_Base(n, magnitude)))
  {
  }

  RealSpace RealSpace::copy() const
  {
    return RealSpace(std::make_shared<RealSpace_Base>(std::static_pointer_cast<RealSpace_Base>(ptr_)->copy()));
  }

}
