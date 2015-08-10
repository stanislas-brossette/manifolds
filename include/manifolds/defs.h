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

#ifndef _MANIFOLDS_DEFS_H_
#define _MANIFOLDS_DEFS_H_

#define EIGEN_RUNTIME_NO_MALLOC

#include <Eigen/Core>

#if defined(_MSC_FULL_VER)
typedef unsigned int uint;
#include <cstdint>
#endif

namespace mnf
{
  typedef Eigen::Ref<Eigen::VectorXd> RefVec;
  typedef Eigen::Ref<const Eigen::VectorXd> ConstRefVec;
  typedef Eigen::Ref<Eigen::MatrixXd> RefMat;
  typedef Eigen::Ref<const Eigen::MatrixXd> ConstRefMat;
  typedef Eigen::Ref<Eigen::VectorXd>::SegmentReturnType Segment;
  typedef Eigen::Ref<const Eigen::VectorXd>::ConstSegmentReturnType ConstSegment;
  typedef Eigen::VectorXd::Index Index;
  const Eigen::IOFormat defaultFormat(4, 0, ", ", "\n", "[", "]");
}

#include <manifolds/manifolds_api.h>

#endif //_MANIFOLDS_DEFS_H_

