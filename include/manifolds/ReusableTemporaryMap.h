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

#ifndef _MANIFOLDS_REUSABLE_TEMPORARY_MAP_H_
#define _MANIFOLDS_REUSABLE_TEMPORARY_MAP_H_

#include <Eigen/Core>

#include <manifolds/defs.h>

namespace mnf
{
  class MANIFOLDS_API ReusableTemporaryMap
  {
  public:
    ReusableTemporaryMap(size_t size = 256);

    //the copy operator only copy the size. The buffer is different
    ReusableTemporaryMap(const ReusableTemporaryMap& other);
    ~ReusableTemporaryMap();

    Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> getMap(Eigen::DenseIndex m, Eigen::DenseIndex n);

  private:
    ReusableTemporaryMap& operator= (const ReusableTemporaryMap&); //We forbid copy

    void reallocate(size_t size);
    void allocate_(size_t size);
    void reallocate_(size_t size);

  private:
    Eigen::aligned_allocator<double> allocator_;
    size_t size_;
    double* buffer_;
  };

  inline Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> ReusableTemporaryMap::getMap(Eigen::DenseIndex m, Eigen::DenseIndex n)
  {
    reallocate(static_cast<size_t>(m*n));
    return Eigen::Map<Eigen::MatrixXd, Eigen::Aligned>(buffer_, m, n);
  }

  inline void ReusableTemporaryMap::reallocate(size_t size)
  {
    if (size > size_)
      reallocate_(size);
  }
}

#endif //_MANIFOLDS_REUSABLE_TEMPORARY_MAP_H_

