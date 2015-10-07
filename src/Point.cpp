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

#include <iostream>

#include <manifolds/defs.h>
#include <manifolds/Point.h>
#include <manifolds/Manifold.h>
#include <manifolds/Manifold_Base.h>
#include <manifolds/mnf_assert.h>

namespace mnf
{

  ConstManifoldOwner::ConstManifoldOwner(const Manifold& m)
    : manifold_(const_cast<Manifold&>(m).shallowCopy())
  {
  }

  const Manifold& ConstManifoldOwner::getManifold() const
  {
    return manifold_;
  }

  //A const cast is necessary to save a Manifold. 
  //We make sure to not modify manifold in this class
ConstSubPoint::ConstSubPoint(const Manifold& m, const ConstRefVec& val)
    : ConstManifoldOwner(m),
      value_(Eigen::Map<Eigen::VectorXd>(const_cast<double*>(val.data()),
                                         val.size())),
      format_(mnf::defaultFormat)
{
  mnf_assert(getManifold().representationDim() == val.size());
}

ConstSubPoint::ConstSubPoint(const ConstSubPoint& other)
    : ConstManifoldOwner(other.getManifold()),
    value_(other.value_)
{
}

ConstSubPoint::~ConstSubPoint() {}

ConstRefVec ConstSubPoint::value() const { return value_; }

const ConstSubPoint ConstSubPoint::operator()(size_t i) const
{
  return ConstSubPoint(getManifold().operator()(i),
                       getManifold().getConstView<R>(value_, i));
}

ConstSegment ConstSubPoint::operator[](size_t i) const
{
  return getManifold().getConstView<R>(value_, i);
}

const Eigen::IOFormat& ConstSubPoint::format() const { return format_; }

std::string ConstSubPoint::toString(std::string& prefix,
                                    const Eigen::IOFormat& fmt) const
{
  return getManifold().toString(value_, prefix, fmt);
}

SubPoint::SubPoint(const Manifold& M, RefVec val) : ConstSubPoint(M, val) {}

SubPoint::SubPoint(const SubPoint& other) : ConstSubPoint(other) {}

RefVec SubPoint::value() { return value_; }

SubPoint SubPoint::operator()(size_t i)
{
  return SubPoint(getManifold().operator()(i), getManifold().getView<R>(value_, i));
}

Segment SubPoint::operator[](size_t i)
{
  return getManifold().getView<R>(value_, i);
}

PointMemory::PointMemory(Index size) : mem_(size) {}

PointMemory::PointMemory(const ConstRefVec& v) : mem_(v) {}

Eigen::VectorXd& PointMemory::getMem()
{
  mnf_assert(mem_.size() > 0);
  return mem_;
}

Point::Point(const Manifold& M)
    : PointMemory(M.representationDim()), SubPoint(M, getMem())
{
}

Point::Point(const Manifold& M, const ConstRefVec& val)
    : PointMemory(val), SubPoint(M, getMem())
{
}

Point::Point(const Point& other)
    : PointMemory(other.value()), SubPoint(other.getManifold(), getMem())
{
}

Point::Point(const ConstSubPoint& other)
    : PointMemory(other.value()), SubPoint(other.getManifold(), getMem())
{
}

Point& Point::increment(const ConstRefVec& v)
{
  getManifold().retractation(value_, value_, v);
  return *this;
}

Point Point::retractation(const ConstRefVec& v) const
{
  Eigen::VectorXd out(getManifold().representationDim());
  getManifold().retractation(out, this->value_, v);
  return getManifold().createPoint(out);
}
void Point::retractation(RefVec out, const ConstRefVec& v) const
{
  getManifold().retractation(out, this->value_, v);
}
void Point::retractation(Point& out, const ConstRefVec& v) const
{
  getManifold().retractation(out.value(), this->value_, v);
}

Eigen::VectorXd Point::pseudoLog(const Point& y) const
{
  Eigen::VectorXd out(getManifold().tangentDim());
  getManifold().pseudoLog(out, this->value_, y.value());
  return out;
}
void Point::pseudoLog(RefVec out, const Point& y) const
{
  getManifold().pseudoLog(out, this->value_, y.value());
}

Eigen::VectorXd Point::pseudoLog0() const
{
  Eigen::VectorXd out(getManifold().tangentDim());
  getManifold().pseudoLog0(out, this->value_);
  return out;
}
void Point::pseudoLog0(RefVec out) const
{
  getManifold().pseudoLog0(out, this->value_);
}

Eigen::VectorXd Point::typicalMagnitude() const
{
  return getManifold().getTypicalMagnitude();
}

Eigen::VectorXd Point::trustMagnitude() const
{
  return getManifold().getTrustMagnitude();
}

Point& Point::operator=(const Point& x)
{
  mnf_assert(getManifold().dim() == x.getManifold().dim());
  mnf_assert(getManifold().representationDim() ==
             x.getManifold().representationDim());
  value_ = x.value();
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

bool Point::isInM(double prec) const { return getManifold().isInM(value_, prec); }
bool Point::isInTxM(const ConstRefVec& v, const double& prec) const
{
  return getManifold().isInTxM(value_, v, prec);
}

Index Point::getDimM() const { return getManifold().dim(); }
Index Point::getTangentDimM() const { return getManifold().tangentDim(); }
Index Point::getRepresentationDimM() const
{
  return getManifold().representationDim();
}

Eigen::MatrixXd Point::diffRetractation() const
{
  return getManifold().diffRetractation(value_);
}

void Point::applyDiffRetractation(RefMat out, const ConstRefMat& in) const
{
  getManifold().applyDiffRetractation(out, in, value_);
}

Eigen::MatrixXd Point::diffPseudoLog0() const
{
  return getManifold().diffPseudoLog0(value_);
}

void Point::applyDiffPseudoLog0(RefMat out, const ConstRefMat& in) const
{
  getManifold().applyDiffPseudoLog0(out, in, value_);
}

void Point::applyTransport(RefMat out, const ConstRefMat& in,
                           const ConstRefVec& v) const
{
  getManifold().applyTransport(out, in, value_, v);
}

void Point::applyInvTransport(RefMat out, const ConstRefMat& in,
                              const ConstRefVec& v) const
{
  getManifold().applyInvTransport(out, in, value_, v);
}

const Point& Point::format(const Eigen::IOFormat& fmt) const
{
  format_ = fmt;
  return *this;
}
}
