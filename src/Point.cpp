#include <iostream>
#include <manifolds/Point.h>
#include <manifolds/Manifold.h>
#include <manifolds/pgs_assert.h>

namespace pgs
{
  ConstSubPoint::ConstSubPoint(const Manifold& M, RefVec val)
    : manifold_(M)
    , value_(val)
  {
    pgs_assert(M.representationDim() == val.size());
    registerPoint();
  }
  
  ConstSubPoint::ConstSubPoint(const Manifold& M, const ConstRefVec& val)
    : manifold_(M)
    , value_(Eigen::Map<Eigen::VectorXd>(const_cast<double*>(val.data()), val.size()))
  {
    pgs_assert(M.representationDim() == val.size());
    registerPoint();
  }

  ConstSubPoint::ConstSubPoint(const ConstSubPoint& other)
    : manifold_(other.manifold_)
    , value_(other.value_)
  {
    registerPoint();
  }

  ConstSubPoint::~ConstSubPoint()
  {
    unregisterPoint();
  }

  ConstRefVec ConstSubPoint::value() const
  {
    return value_;
  }

  const ConstSubPoint ConstSubPoint::operator()(size_t i) const
  {
    return ConstSubPoint(manifold_(i), manifold_.getConstView<R>(value_, i));
  }

  ConstSegment ConstSubPoint::operator[](size_t i) const
  {
    return manifold_.getConstView<R>(value_, i);
  }

  const Manifold& ConstSubPoint::getManifold() const
  {
    return manifold_;
  }

  std::string ConstSubPoint::toString(std::string& prefix) const
  {
    return manifold_.toString(value_, prefix);
  }

  void ConstSubPoint::registerPoint()
  {
    this->manifold_.incrementRefCounter();
  }

  void ConstSubPoint::unregisterPoint()
  {
    this->manifold_.decrementRefCounter();
  }



  SubPoint::SubPoint(const Manifold& M, RefVec val)
    : ConstSubPoint(M, val)
  {
  }

  SubPoint::SubPoint(const SubPoint& other)
    : ConstSubPoint(other)
  {
  }

  RefVec SubPoint::value()
  {
    return value_;
  }

  SubPoint SubPoint::operator()(size_t i)
  {
    return SubPoint(manifold_(i), manifold_.getView<R>(value_, i));
  }

  Segment SubPoint::operator[](size_t i)
  {
    return manifold_.getView<R>(value_, i);
  }



  PointMemory::PointMemory(Index size) : mem_(size) {}

  PointMemory::PointMemory(const ConstRefVec& v) : mem_(v) {}

  Eigen::VectorXd& PointMemory::getMem()
  {
    pgs_assert(mem_.size() > 0);
    return mem_;
  }



  Point::Point(const Manifold& M)
    : PointMemory(M.representationDim())
    , SubPoint(M, getMem())
  {
  }

  Point::Point(const Manifold& M, const ConstRefVec& val)
    : PointMemory(val)
    , SubPoint(M, getMem())
  {
  }

  Point::Point(const Point& other)
    : PointMemory(other.value())
    , SubPoint(other.getManifold(), getMem())
  {
  }

  Point::Point(const ConstSubPoint& other)
    : PointMemory(other.value())
    , SubPoint(other.getManifold(), getMem())
  {
  }

  Point& Point::increment(const ConstRefVec& v)
  {
    manifold_.retractation(value_, value_, v);
    return *this;
  }

  Point Point::retractation(const ConstRefVec& v)
  {
    Eigen::VectorXd out(manifold_.representationDim());
    manifold_.retractation(out, this->value_, v);
    return manifold_.createPoint(out);
  }
   
  RefVec Point::pseudoLog(const Point& y)
  {
    Eigen::VectorXd out(manifold_.tangentDim());
    manifold_.pseudoLog(out, this->value_, y.value());
    return out;
  }
  
  RefVec Point::pseudoLog0()
  {
    Eigen::VectorXd out(manifold_.tangentDim());
    manifold_.pseudoLog0(out, this->value_);
    return out;
  }

  Point & Point::operator=(const Point& x)
  {
    pgs_assert(this->manifold_.dim() == x.getManifold().dim());
    pgs_assert(this->manifold_.representationDim() == x.getManifold().representationDim());
    this->value_ = x.value();
    return *this;
  }
  
  Point operator+(const Point& x, const ConstRefVec& v)
  {
    return x.getManifold().createPoint(x.value()).increment(v);
  }

  Eigen::VectorXd operator-(const Point& x, const Point& y)
  {
    Eigen::VectorXd output(x.getManifold().dim());
    x.getManifold().pseudoLog(output, y.value(), x.value());
    return output;
  }



}
