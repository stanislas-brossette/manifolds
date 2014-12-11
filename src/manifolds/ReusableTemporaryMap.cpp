#include <pgsolver/manifolds/ReusableTemporaryMap.h>

namespace pgs
{
  ReusableTemporaryMap::ReusableTemporaryMap(size_t size)
    : size_(0)
    , buffer_(0x0)
  {
    assert(size > 0 && "size must be at least one");
    allocate_(size);
  }

  ReusableTemporaryMap::ReusableTemporaryMap(const ReusableTemporaryMap& other)
    : size_(0)
    , buffer_(0x0)
  {
    allocate_(other.size_);
  }

  ReusableTemporaryMap::~ReusableTemporaryMap()
  {
    assert(buffer_ != 0x0);
    allocator_.deallocate(buffer_, size_);
  }

  void ReusableTemporaryMap::allocate_(size_t size)
  {
    assert(buffer_ == 0x0);
    buffer_ = allocator_.allocate(size);
    size_ = size;
  }

  void ReusableTemporaryMap::reallocate_(size_t size)
  {
    size_t newSize = size_;
    while (newSize < size)
      newSize *= 2;
    allocator_.deallocate(buffer_, size_);
    buffer_ = allocator_.allocate(newSize);
    size_ = newSize;
  }
}
