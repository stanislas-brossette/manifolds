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
  Manifold::Manifold(Index dimension, Index tangentDimension, Index representationDimension)
    :   
      manifoldBase_(make_shared<Manifold_Base>(dimension, tangentDimension, representationDimension))
  {}

  Manifold(std::shared_ptr<Manifold_Base> m)
    :
      manifoldBase_(m)
  {}

  Manifold::~Manifold()
  {
  }

  Point Manifold::createPoint() const
  {
    return manifoldBase_->createPoint();
  }

  Point Manifold::createPoint(const ConstRefVec& val) const
  {
    return manifoldBase_->createPoint(val);
  }

  bool Manifold::isInM(const Eigen::VectorXd& val, double prec) const
  {
    return manifoldBase_->isInM(val, prec);
  }

  void Manifold::forceOnM(RefVec out, const ConstRefVec& in) const
  {
    return manifoldBase_->forceOnM(out, in);
  }

  Point Manifold::getZero() const
  {
    return manifoldBase_->getZero();
  }

  Point Manifold::createRandomPoint(double coeff) const
  {
    return manifoldBase_->createRandomPoint(coeff);
  }

  void Manifold::createRandomPoint(RefVec out, double coeff) const
  {
    return manifoldBase_->createRandomPoint(out, coeff);
  }

  Index Manifold::dim() const
  {
    return manifoldBase_->dim();
  }

  Index Manifold::tangentDim() const
  {
    return manifoldBase_->tangentDim();
  }

  Index Manifold::representationDim() const
  {
    return manifoldBase_->representationDim();
  }

  void Manifold::display(const std::string& prefix) const
  {
    return manifoldBase_->display(prefix);
  }

  void Manifold::retractation(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    return manifoldBase_->retractation(out, x, v);
  }

  void Manifold::pseudoLog(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    return manifoldBase_->pseudoLog(out, x, y);
  }

  void Manifold::pseudoLog0(RefVec out, const ConstRefVec& x) const
  {
    return manifoldBase_->pseudoLog0(out, x);
  }

  Eigen::MatrixXd Manifold::diffRetractation(const ConstRefVec& x) const
  {
    return manifoldBase_->diffRetractation(x);
  }

  void Manifold::applyDiffRetractation( RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    return manifoldBase_->applyDiffRetractation(out, in, x);
  }

  Eigen::MatrixXd Manifold::diffPseudoLog0(const ConstRefVec& x) const
  {
    return manifoldBase_->diffPseudoLog0(x);
  }

  void Manifold::applyDiffPseudoLog0(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    return manifoldBase_->applyDiffPseudoLog0(out, in, x);
  }

  void Manifold::applyTransport(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    return manifoldBase_->applyTransport(out, in, x, v);
  }

  void Manifold::applyInvTransport(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    return manifoldBase_->applyInvTransport(out, in, x, v);
  }

  void Manifold::applyInvTransportOnTheRight(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    return manifoldBase_->applyInvTransportOnTheRight(out, in, x, v);
  }

  void Manifold::setZero(RefVec out) const
  {
    return manifoldBase_->setZero(out);
  }

  void Manifold::tangentConstraint(RefMat out, const ConstRefVec& x) const
  {
    return manifoldBase_->tangentConstraint(out, x);
  }

  bool Manifold::isInTxM(const ConstRefVec& x, const ConstRefVec& v, const double& prec) const
  {
    return manifoldBase_->isInTxM(x, v, prec);
  }

  void Manifold::forceOnTxM(RefVec out, const ConstRefVec& in, const ConstRefVec&x) const
  {
    return manifoldBase_->forceOnTxM(out, in, x);
  }

  void Manifold::getIdentityOnTxM(RefMat out, const ConstRefVec& x) const
  {
    return manifoldBase_->getIdentityOnTxM(out, x);
  }

  void Manifold::limitMap(RefVec out) const
  {
    return manifoldBase_->limitMap(out);
  }

  void Manifold::getTypicalMagnitude(RefVec out) const
  {
    return manifoldBase_->getTypicalMagnitude(out);
  }

  Eigen::VectorXd Manifold::getTypicalMagnitude() const
  {
    return manifoldBase_->VectorXd Manifold::getTypicalMagnitude();
  }

  void Manifold::getTrustMagnitude(RefVec out) const
  {
    return manifoldBase_->getTrustMagnitude(out);
  }

  Eigen::VectorXd Manifold::getTrustMagnitude() const
  {
    return manifoldBase_->VectorXd Manifold::getTrustMagnitude();
  }

  void Manifold::lock() const
  {
    return manifoldBase_->lock();
  }

  bool Manifold::isLocked() const
  {
    return manifoldBase_->isLocked();
  }

  const std::string& Manifold::name() const 
  {
    return manifoldBase_->name();
  }

  std::string& Manifold::name() 
  {
    return manifoldBase_->name();
  }

  long Manifold::getInstanceId() const
  {
    return manifoldBase_->getInstanceId();
  }

  bool Manifold::isSameType(const Manifold& other) const
  {
    return manifoldBase_->isSameType(other);
  }
}
