#include "Point.h"
#include "Manifold.h"
#include <assert.h>

namespace pgs
{
  Point::Point(const Manifold& M)
    : manifold_(M)
    , value_(M.representationDim())
  {
  }


  Point::Point(const Manifold& M, const Eigen::VectorXd& val)
    : manifold_(M)
    , value_(val)
  {
    assert(static_cast<int> (M.representationDim()) == val.size());
  }

  
  Point& Point::increment(const Eigen::VectorXd& v)
  {
    manifold_.plus(value_, value_, v);
    return *this;
  }


  Point Point::operator()(size_t i) const
  {
    return Point(manifold_(i), manifold_.getValueConst(value_, i));
  }

  const Eigen::VectorXd& Point::value() const
  {
    return value_;
  }


  ConstSegment Point::operator[](size_t i) const
  {
    return manifold_.getValueConst(value_, i);
  }

  Segment Point::operator[](size_t i)
  {
    return manifold_.getValue(value_, i);
  }
  
  const Manifold& Point::getManifold() const
  {
    return manifold_;
  }

  Point operator+(const Point& x, const Eigen::Ref<const Eigen::VectorXd>& v)
  {
    return x.getManifold().createPoint(x.value()).increment(v);
  }
}
