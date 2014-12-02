#include <pgsolver/Point.h>
#include <pgsolver/Manifold.h>
#include <assert.h>
#include <iostream>

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
    assert(M.representationDim() == val.size());
  }

  Point::Point(const Point& other)
    : manifold_(other.manifold_)
    , value_(other.value_)
  {
    other.manifold_.incrementRefCounter();
  }

  Point::~Point()
  {
    manifold_.decrementRefCounter();
  }

  
  Point& Point::increment(const Eigen::VectorXd& v)
  {
    manifold_.plus(value_, value_, v);
    return *this;
  }
  
  Eigen::VectorXd Point::invMap() const
  {
    Eigen::VectorXd res(manifold_.dim());
    manifold_.minus(res,value_,manifold_.getIdentity().value());
    return res;
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

  Point operator+(const Point& x, ConstRefVec& v)
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
    assert(this->manifold_.dim() == x.manifold_.dim());
    assert(this->manifold_.representationDim() == x.manifold_.representationDim());
    this->value_ = x.value();
    return *this;
  }

  std::string Point::toString(std::string& prefix) const
  {
    return manifold_.toString(value_, prefix);
  }
}
