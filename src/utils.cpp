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
#include <Eigen/Geometry>

namespace mnf
{

namespace utils
{
  void hat(Eigen::Matrix3d& M, const Eigen::Vector3d& v)
  {
    M << 0    , -v(2), v(1) ,
         v(2) , 0    , -v(0),
         -v(1), v(0) , 0    ;
  }
  void hat2(Eigen::Matrix3d& M, const Eigen::Vector3d& v)
  {
    double v0v1 = v(0)*v(1);
    double v1v2 = v(1)*v(2);
    double v2v0 = v(2)*v(0);
    double v0_2 = v(0)*v(0);
    double v1_2 = v(1)*v(1);
    double v2_2 = v(2)*v(2);

    M << -v1_2-v2_2, v0v1, v2v0,
         v0v1, -v0_2-v2_2, v1v2,
         v2v0, v1v2, -v0_2-v1_2;
  }

  Eigen::Matrix3d computeRotBetweenVec(const Eigen::Vector3d& x, const Eigen::Vector3d& y)
  {
    Eigen::Vector3d acrossb = x.cross(y);
    double adotb = x.dot(y);
    assert(1+adotb > 1e-12);
    Eigen::Matrix3d acrossbHat;
    Eigen::Matrix3d acrossbHat2;
    utils::hat(acrossbHat, acrossb);
    utils::hat2(acrossbHat2, acrossb);
    Eigen::Matrix3d R;
    R = Eigen::Matrix3d::Identity() + acrossbHat + acrossbHat2/(1+adotb);
    return R;
  }

  bool set_is_malloc_allowed(bool allow __attribute__((unused)) )
  {
#ifdef EIGEN_RUNTIME_NO_MALLOC
    return Eigen::internal::set_is_malloc_allowed(allow);
#else
    eigen_assert(false && "you can't call this function if Manifold was compiled without the flag EIGEN_RUNTIME_NO_MALLOC");
    return true;
#endif
  }
}

}
