#include <stdexcept>
#include "Manifold.h"

namespace pgs
{
  Manifold::Manifold(size_t dimension, size_t representationDimension)
    : dimension_(dimension)
    , representationDim_(representationDimension)
    , lock_(false)
  {
  }

  Point Manifold::createPoint() const
  {
    lock();
    return Point(*this);
  }

  Point Manifold::createPoint(const Eigen::VectorXd& val) const
  {
    lock();
    return Point(*this, val);
  }

  Point Manifold::getIdentity() const
  {
    lock();
    Eigen::VectorXd id(representationDim_);
    setIdentity(id);
    return Point(*this, id);
  }

  size_t Manifold::dim() const
  {
    return dimension_;
  }

  size_t Manifold::representationDim() const
  {
    return representationDim_;
  }

  void Manifold::plus(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v) const
  {
    assert(out.size() == static_cast<int> (representationDim_));
    assert(x.size() == static_cast<int> (representationDim_));
    assert(v.size() == static_cast<int> (dimension_));
    plus_(out, x, v);
  }


  void Manifold::setDimension(size_t d)
  {
    testLock();
    dimension_ = d;
  }

  void Manifold::setRepresentationDimension(size_t rd)
  {
    testLock();
    representationDim_ = rd;
  }


  void Manifold::setIdentity(Eigen::Ref<Eigen::VectorXd> out) const
  {
    assert(out.size() == static_cast<int> (representationDim_));
    setIdentity_(out);
  }
  

  void Manifold::lock() const
  {
    lock_ = true;
  }

  void Manifold::testLock() const
  {
    if (lock_)
      throw std::runtime_error("Either a point or a compound manifold is relying on this manifold, you can't modify it anymore.");
  }
}
