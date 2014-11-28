#include <stdexcept>
#include "Manifold.h"

namespace pgs
{
  Manifold::Manifold(Index dimension, Index representationDimension)
    : dimension_(dimension)
    , representationDim_(representationDimension)
    , lock_(false)
  {
    assert(dimension>=0 && "Negative dimension not accepted");
    assert(representationDimension>=0 && "Negative dimension not accepted");
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

  Index Manifold::dim() const
  {
    return dimension_;
  }

  Index Manifold::representationDim() const
  {
    return representationDim_;
  }

  void Manifold::plus(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v) const
  {
    assert(out.size() ==  representationDim_);
    assert(x.size() == representationDim_);
    assert(v.size() == dimension_);
    plus_(out, x, v);
  }


  void Manifold::setDimension(Index d)
  {
    assert(d>0 && "Negative dimension not accepted");
    testLock();
    dimension_ = d;
  }

  void Manifold::setRepresentationDimension(Index rd)
  {
    assert(rd>0 && "Negative dimension not accepted");
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
