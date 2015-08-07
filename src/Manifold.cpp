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

#include <stdexcept>
#include <manifolds/Manifold.h>
#include <manifolds/mnf_assert.h>

namespace mnf
{
  long Manifold::manifoldCounter_ = 0;

  Manifold::Manifold(Index dimension, Index tangentDimension, Index representationDimension)
    : dimension_(dimension)
    , tangentDim_(tangentDimension)
    , representationDim_(representationDimension)
    , lock_(false)
  {
    mnf_assert(0 <= dimension && "Negative dimension not accepted");
    mnf_assert(dimension <= tangentDimension);
    mnf_assert(tangentDimension <= representationDimension);
    this->instanceId_ = manifoldCounter_++;
  }

  Manifold::~Manifold()
  {
  }

  Point Manifold::createPoint() const
  {
    mnf_assert(isValid() || seeMessageAbove());
    lock();
    return Point(*this);
  }

  Point Manifold::createPoint(const ConstRefVec& val) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    if (isInM(val))
    {
      lock();
      return Point(*this, val);
    }
    else
    {
      throw std::runtime_error("Bad Point Initialization");
    }
  }

  bool Manifold::isInM(const Eigen::VectorXd& val, const double& prec) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(val.size() == representationDim());
    return isInM_(val, prec);
  }

  void Manifold::forceOnM(RefVec out, const ConstRefVec& in) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(out.size() == representationDim());
    mnf_assert(in.size() == representationDim());
    return forceOnM_(out, in);
  }


  Point Manifold::getZero() const
  {
    mnf_assert(isValid() || seeMessageAbove());
    lock();
    Eigen::VectorXd id(representationDim_);
    setZero(id);
    return Point(*this, id);
  }

  Point Manifold::createRandomPoint(double coeff) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    lock();
    Eigen::VectorXd val(representationDim_);
    createRandomPoint(val, coeff);
    return Point(*this, val);
  }

  void Manifold::createRandomPoint(RefVec out, double coeff) const
  {
    mnf_assert(out.size() == representationDim_ && "wrong dimension in Manifold::createRandomPoint");
    createRandomPoint_(out, coeff);
  }

  Index Manifold::dim() const
  {
    mnf_assert(isValid() || seeMessageAbove());
    return dimension_;
  }

  Index Manifold::tangentDim() const
  {
    mnf_assert(isValid() || seeMessageAbove());
    return tangentDim_;
  }

  Index Manifold::representationDim() const
  {
    mnf_assert(isValid() || seeMessageAbove());
    return representationDim_;
  }

  void Manifold::display(std::string prefix) const
  {
    std::cout << prefix << name() << std::endl;
  }

  void Manifold::retractation(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(out.size() == representationDim_);
    mnf_assert(x.size() == representationDim_);
    mnf_assert(v.size() == tangentDim_);
    mnf_assert(isInTxM(x, v) && "Wrong tangent vector provided to retractation");
    retractation_(out, x, v);
  }
  //void Manifold::retractation(RefVec out, const Point& x, const ConstRefVec& v) const
  //{
  //  retractation( out, x.value(), v);
  //}
  //Point Manifold::retractation(const ConstRefVec& x, const ConstRefVec& v) const
  //{
  //  Eigen::VectorXd out(representationDim_);
  //  retractation(out, x, v);
  //  return createPoint(out);
  //}
  //Point Manifold::retractation(const Point& x, const ConstRefVec& v) const
  //{
  //  return retractation(x.value(), v);
  //}

  void Manifold::pseudoLog(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(out.size() == tangentDim_);
    mnf_assert(x.size() == representationDim_);
    mnf_assert(y.size() == representationDim_);
    pseudoLog_(out, x, y);
  }
  //void Manifold::pseudoLog(RefVec out, const Point& x, const Point& y) const
  //{
  //  pseudoLog(out, x.value(), y.value());
  //}
  //RefVec Manifold::pseudoLog(const ConstRefVec& x, const ConstRefVec& y) const
  //{
  //  Eigen::VectorXd out(tangentDim_);
  //  pseudoLog(out, x, y);
  //  return out;
  //}
  //RefVec Manifold::pseudoLog(const Point& x, const Point& y) const
  //{
  //  Eigen::VectorXd out(tangentDim_);
  //  pseudoLog(out, x.value(), y.value());
  //  return out;
  //}

  void Manifold::pseudoLog0(RefVec out, const ConstRefVec& x) const
  {
    mnf_assert(out.size() == tangentDim_);
    mnf_assert(x.size() == representationDim_);
    pseudoLog0_(out, x);
  }
  //void Manifold::pseudoLog0(RefVec out, const Point& x) const
  //{
  //  pseudoLog0(out, x.value());
  //}
  //RefVec Manifold::pseudoLog0(const ConstRefVec& x) const
  //{
  //  Eigen::VectorXd out(tangentDim_);
  //  pseudoLog0(out, x);
  //  return out;
  //}
  //RefVec Manifold::pseudoLog0(const Point& x) const
  //{
  //  Eigen::VectorXd out(tangentDim_);
  //  pseudoLog0(out, x.value());
  //  return out;
  //}

  Eigen::MatrixXd Manifold::diffRetractation(const ConstRefVec& x) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(x.size() == representationDim_);
    return diffRetractation_(x);
  }

  void Manifold::applyDiffRetractation(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(in.cols() == representationDim_);
    mnf_assert(out.cols() == tangentDim_);
    mnf_assert(in.rows() == out.rows());
    mnf_assert(x.size() == representationDim_);
    applyDiffRetractation_(out, in, x);
  }

  Eigen::MatrixXd Manifold::diffPseudoLog0(const ConstRefVec& x) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(x.size() == representationDim_);
    return diffPseudoLog0_(x);
  }

  void Manifold::applyDiffPseudoLog0(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(out.cols() == representationDim_);
    mnf_assert(in.cols() == tangentDim_);
    mnf_assert(in.rows() == out.rows());
    mnf_assert(x.size() == representationDim_);
    applyDiffPseudoLog0_(out, in, x);
  }

  void Manifold::applyTransport(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(in.rows() == tangentDim_);
    mnf_assert(out.rows() == tangentDim_);
    mnf_assert(in.cols() == out.cols());
    mnf_assert(x.size() == representationDim());
    mnf_assert(v.size() == tangentDim_);
    mnf_assert(isInTxM(x, v));

    applyTransport_(out, in, x, v);
  }

  void Manifold::applyInvTransport(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(in.rows() == tangentDim_);
    mnf_assert(out.rows() == tangentDim_);
    mnf_assert(in.cols() == out.cols());
    mnf_assert(x.size() == representationDim());
    mnf_assert(v.size() == tangentDim_);
    mnf_assert(isInTxM(x, v));

    applyInvTransport_(out, in, x, v);
  }

  void Manifold::applyInvTransportOnTheRight(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(in.cols() == tangentDim_);
    mnf_assert(out.cols() == tangentDim_);
    mnf_assert(in.rows() == out.rows());
    mnf_assert(x.size() == representationDim());
    mnf_assert(v.size() == tangentDim_);
    mnf_assert(isInTxM(x, v));

    applyInvTransportOnTheRight_(out, in, x, v);
  }

  void Manifold::setDimension(Index d)
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(d>=0 && "Negative dimension not accepted");
    testLock();
    dimension_ = d;
  }

  void Manifold::setTangentDimension(Index td)
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(td >= 0 && "Negative dimension not accepted");
    testLock();
    tangentDim_ = td;
  }

  void Manifold::setRepresentationDimension(Index rd)
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(rd>=0 && "Negative dimension not accepted");
    testLock();
    representationDim_ = rd;
  }


  void Manifold::setZero(RefVec out) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(out.size() == representationDim_);
    setZero_(out);
  }

  void Manifold::tangentConstraint(RefMat out, const ConstRefVec& x) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(out.rows() == tangentDim_ - dimension_);
    mnf_assert(out.cols() == tangentDim_);
    mnf_assert(x.size() == representationDim());
    tangentConstraint_(out, x);
  }

  bool Manifold::isInTxM(const ConstRefVec& x, const ConstRefVec& v, const double& prec) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(v.size() == tangentDim_);
    mnf_assert(x.size() == representationDim());
    return isInTxM_(x, v, prec);
  }

  void Manifold::forceOnTxM(RefVec out, const ConstRefVec& in, const ConstRefVec&x) const
  {
    mnf_assert(isValid() || seeMessageAbove());
    mnf_assert(out.size() == tangentDim_);
    mnf_assert(x.size() == representationDim());
    mnf_assert(in.size() == tangentDim_);
    forceOnTxM_(out, in, x);
  }

  void Manifold::getIdentityOnTxM(RefMat out, const ConstRefVec& x) const
  {
    mnf_assert(out.rows() == tangentDim_);
    mnf_assert(out.cols() == tangentDim_);
    mnf_assert(x.size() == representationDim());
    getIdentityOnTxM_(out, x);
  }

  void Manifold::limitMap(RefVec out) const
  {
    mnf_assert(out.size() == tangentDim_);
    limitMap_(out);
  }

  void Manifold::getTypicalMagnitude(RefVec out) const
  {
    mnf_assert(out.size() == tangentDim_);
    getTypicalMagnitude_(out);
  }

  Eigen::VectorXd Manifold::getTypicalMagnitude() const
  {
    Eigen::VectorXd out(tangentDim_);
    getTypicalMagnitude_(out);
    return out;
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

  bool Manifold::isLocked() const
  {
    return lock_;
  }

  const std::string& Manifold::name() const { return name_;}

  std::string& Manifold::name() { return name_;}

  long Manifold::getInstanceId() const
  {
    return this->instanceId_;
  }

  bool Manifold::isSameType(const Manifold& other) const
  {
    return this->getTypeId() == other.getTypeId();
  }

  Manifold* Manifold::copyManifold(const Manifold& m)
  {
    Manifold* copy = m.getNewCopy();

    copy->instanceId_ = m.instanceId_;

    return copy;
  }
}
