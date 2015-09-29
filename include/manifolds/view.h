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

#include <Eigen/Core>
#include <manifolds/defs.h>

namespace mnf
{
  enum eDimension
  {
    R,  //representation space
    T,  //tangent space
    F  //full space
  };

  template<int Dr, int Dc> struct ViewReturnType { typedef Eigen::Block<RefMat> Type; };
  template<int Dr, int Dc> struct ConstViewReturnType { typedef const Eigen::Block<ConstRefMat> Type; };
  template<int Dc> struct ViewReturnType<F, Dc> { typedef RefMat::ColsBlockXpr Type; };
  template<int Dc> struct ConstViewReturnType<F, Dc> { typedef ConstRefMat::ConstColsBlockXpr Type; };
  template<int Dr> struct ViewReturnType<Dr, F> { typedef RefMat::RowsBlockXpr Type; };
  template<int Dr> struct ConstViewReturnType<Dr, F> { typedef ConstRefMat::ConstRowsBlockXpr Type; };
  template<> struct ViewReturnType<F, F> { typedef RefMat Type; };
  template<> struct ConstViewReturnType<F, F> { typedef ConstRefMat Type; };
}

