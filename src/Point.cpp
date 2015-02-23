#include <iostream>
#include <manifolds/Point.h>
#include <manifolds/Manifold.h>
#include <manifolds/pgs_assert.h>

namespace pgs
{
  Point::Point(const Manifold& M)
    : manifold_(M)
    , value_(M.representationDim())
  {
    registerPoint();
  }

  Point::Point(const Manifold& M, const Eigen::VectorXd& val)
    : manifold_(M)
    , value_(val)
  {
    pgs_assert(M.representationDim() == val.size());
    registerPoint();
  }

  Point::Point(const Point& other)
    : manifold_(other.manifold_)
    , value_(other.value_)
  {
    registerPoint();
  }

  Point::~Point()
  {
    unregisterPoint();
  }


  Point& Point::increment(const Eigen::VectorXd& v)
  {
    manifold_.plus(value_, value_, v);
    return *this;
  }

  Point Point::operator()(size_t i) const
  {
    return Point(manifold_(i), manifold_.getConstView<R>(value_, i));
  }

  const Eigen::VectorXd& Point::value() const
  {
    return value_;
  }

  ConstSegment Point::operator[](size_t i) const
  {
    pgs_assert(manifold_.isValid() || manifold_.seeMessageAbove());
    return manifold_.getConstView<R>(value_, i);
  }

  Segment Point::operator[](size_t i)
  {
    return manifold_.getView<R>(value_, i);
  }

  const Manifold& Point::getManifold() const
  {
    return manifold_;
  }

  Point operator+(const Point& x, const ConstRefVec& v)
  {
    return x.getManifold().createPoint(x.value()).increment(v);
  }

  Eigen::VectorXd operator-(const Point& x, const Point& y)
  {
    Eigen::VectorXd output(x.getManifold().dim());
    x.getManifold().minus(output, x.value(), y.value());
    return output;
  }

  Point & Point::operator=(const Point& x)
  {
    pgs_assert(this->manifold_.dim() == x.manifold_.dim());
    pgs_assert(this->manifold_.representationDim() == x.manifold_.representationDim());
    this->value_ = x.value();
    return *this;
  }

  std::string Point::toString(std::string& prefix) const
  {
    return manifold_.toString(value_, prefix);
  }

  void Point::registerPoint()
  {
    this->manifold_.incrementRefCounter();
  }

  void Point::unregisterPoint()
  {
    this->manifold_.decrementRefCounter();
  }
}
