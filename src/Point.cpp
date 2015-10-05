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
ConstSubPoint::ConstSubPoint(ConstManifold M, const ConstRefVec& val)
    : manifold_(M),
      value_(Eigen::Map<Eigen::VectorXd>(const_cast<double*>(val.data()),
                                         val.size())),
      format_(mnf::defaultFormat)
{
  mnf_assert(manifold_.representationDim() == val.size());
}

ConstSubPoint::ConstSubPoint(const ConstSubPoint& other)
    : manifold_(other.manifold_), value_(other.value_)
{
}

ConstSubPoint::~ConstSubPoint() {}

ConstRefVec ConstSubPoint::value() const { return value_; }

const ConstSubPoint ConstSubPoint::operator()(size_t i) const
{
  return ConstSubPoint(manifold_.operator()(i),
                       manifold_.getConstView<R>(value_, i));
}

ConstSegment ConstSubPoint::operator[](size_t i) const
{
  return manifold_.getConstView<R>(value_, i);
}

ConstManifold ConstSubPoint::getManifold() const
{
  return ConstManifold(manifold_);
}

const Eigen::IOFormat& ConstSubPoint::format() const { return format_; }

std::string ConstSubPoint::toString(std::string& prefix,
                                    const Eigen::IOFormat& fmt) const
{
  return manifold_.toString(value_, prefix, fmt);
}

SubPoint::SubPoint(ConstManifold M, RefVec val) : ConstSubPoint(M, val) {}

SubPoint::SubPoint(const SubPoint& other) : ConstSubPoint(other) {}

RefVec SubPoint::value() { return value_; }

SubPoint SubPoint::operator()(size_t i)
{
  return SubPoint(manifold_.operator()(i), manifold_.getView<R>(value_, i));
}

Segment SubPoint::operator[](size_t i)
{
  return manifold_.getView<R>(value_, i);
}

PointMemory::PointMemory(Index size) : mem_(size) {}

PointMemory::PointMemory(const ConstRefVec& v) : mem_(v) {}

Eigen::VectorXd& PointMemory::getMem()
{
  mnf_assert(mem_.size() > 0);
  return mem_;
}

Point::Point(ConstManifold M)
    : PointMemory(M.representationDim()), SubPoint(M, getMem())
{
}

Point::Point(ConstManifold M, const ConstRefVec& val)
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
  manifold_.retractation(value_, value_, v);
  return *this;
}

Point Point::retractation(const ConstRefVec& v) const
{
  Eigen::VectorXd out(manifold_.representationDim());
  manifold_.retractation(out, this->value_, v);
  return manifold_.createPoint(out);
}
void Point::retractation(RefVec out, const ConstRefVec& v) const
{
  manifold_.retractation(out, this->value_, v);
}
void Point::retractation(Point& out, const ConstRefVec& v) const
{
  manifold_.retractation(out.value(), this->value_, v);
}

Eigen::VectorXd Point::pseudoLog(const Point& y) const
{
  Eigen::VectorXd out(manifold_.tangentDim());
  manifold_.pseudoLog(out, this->value_, y.value());
  return out;
}
void Point::pseudoLog(RefVec out, const Point& y) const
{
  manifold_.pseudoLog(out, this->value_, y.value());
}

Eigen::VectorXd Point::pseudoLog0() const
{
  Eigen::VectorXd out(manifold_.tangentDim());
  manifold_.pseudoLog0(out, this->value_);
  return out;
}
void Point::pseudoLog0(RefVec out) const
{
  manifold_.pseudoLog0(out, this->value_);
}

Eigen::VectorXd Point::typicalMagnitude() const
{
  return manifold_.getTypicalMagnitude();
}

Eigen::VectorXd Point::trustMagnitude() const
{
  return manifold_.getTrustMagnitude();
}

Point& Point::operator=(const Point& x)
{
  mnf_assert(this->manifold_.dim() == x.getManifold().dim());
  mnf_assert(this->manifold_.representationDim() ==
             x.getManifold().representationDim());
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

bool Point::isInM(double prec) const { return manifold_.isInM(value_, prec); }
bool Point::isInTxM(const ConstRefVec& v, const double& prec) const
{
  return manifold_.isInTxM(value_, v, prec);
}

Index Point::getDimM() const { return manifold_.dim(); }
Index Point::getTangentDimM() const { return manifold_.tangentDim(); }
Index Point::getRepresentationDimM() const
{
  return manifold_.representationDim();
}

Eigen::MatrixXd Point::diffRetractation() const
{
  return manifold_.diffRetractation(value_);
}

void Point::applyDiffRetractation(RefMat out, const ConstRefMat& in) const
{
  manifold_.applyDiffRetractation(out, in, value_);
}

Eigen::MatrixXd Point::diffPseudoLog0() const
{
  return manifold_.diffPseudoLog0(value_);
}

void Point::applyDiffPseudoLog0(RefMat out, const ConstRefMat& in) const
{
  manifold_.applyDiffPseudoLog0(out, in, value_);
}

void Point::applyTransport(RefMat out, const ConstRefMat& in,
                           const ConstRefVec& v) const
{
  manifold_.applyTransport(out, in, value_, v);
}

void Point::applyInvTransport(RefMat out, const ConstRefMat& in,
                              const ConstRefVec& v) const
{
  manifold_.applyInvTransport(out, in, value_, v);
}

const Point& Point::format(const Eigen::IOFormat& fmt) const
{
  format_ = fmt;
  return *this;
}
}
