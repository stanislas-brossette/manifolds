#ifndef _MANIFOLDS_REUSABLE_TEMPORARY_MAP_H_
#define _MANIFOLDS_REUSABLE_TEMPORARY_MAP_H_

#include <Eigen/Core>

#include <manifolds/defs.h>
#include <pgs_api.h>

namespace pgs
{
  class PGS_API ReusableTemporaryMap
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

