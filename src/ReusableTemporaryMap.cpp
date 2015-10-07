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

#include <manifolds/ReusableTemporaryMap.h>
#include <manifolds/mnf_assert.h>

namespace mnf
{
namespace internal
{
ReusableTemporaryMap::ReusableTemporaryMap(size_t size) : size_(0), buffer_(0x0)
{
  mnf_assert(size > 0 && "size must be at least one");
  allocate_(size);
}

ReusableTemporaryMap::ReusableTemporaryMap(const ReusableTemporaryMap& other)
    : size_(0), buffer_(0x0)
{
  allocate_(other.size_);
}

ReusableTemporaryMap::~ReusableTemporaryMap()
{
  mnf_assert(buffer_ != 0x0);
  allocator_.deallocate(buffer_, size_);
}

void ReusableTemporaryMap::allocate_(size_t size)
{
  mnf_assert(buffer_ == 0x0);
  buffer_ = allocator_.allocate(size);
  size_ = size;
}

void ReusableTemporaryMap::reallocate_(size_t size)
{
  size_t newSize = size_;
  while (newSize < size) newSize *= 2;
  allocator_.deallocate(buffer_, size_);
  buffer_ = allocator_.allocate(newSize);
  size_ = newSize;
}
}
}
