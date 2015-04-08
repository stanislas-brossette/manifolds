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
#include <manifolds/utils.h>

#include <Eigen/Core>

namespace utils
{
  bool set_is_malloc_allowed(bool allow)
  {
#ifdef EIGEN_RUNTIME_NO_MALLOC
    return Eigen::internal::set_is_malloc_allowed(allow);
#else
    eigen_assert(false && "you can't call this function if EigenQP was compiled without the flag EIGEN_RUNTIME_NO_MALLOC");
    return true;
#endif
  }
}
