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

#include <manifolds/Point.h>
#include <manifolds/Manifold.h>
#include <manifolds/mnf_assert.h>

namespace mnf
{
ConstManifold::ConstManifold(std::shared_ptr<const Manifold_Base> m)
    : ptr_(std::const_pointer_cast<Manifold_Base>(m))
{
}

ConstManifold::ConstManifold(std::shared_ptr<Manifold_Base> m) : ptr_(m) {}

ConstManifold::~ConstManifold() {}

std::shared_ptr<const Manifold_Base> ConstManifold::ptr() const
{
  return std::const_pointer_cast<const Manifold_Base>(ptr_);
}

Point ConstManifold::createPoint() const { return ptr_->createPoint(); }

Point ConstManifold::createPoint(const ConstRefVec& val) const
{
  return ptr_->createPoint(val);
}

bool ConstManifold::isInM(const Eigen::VectorXd& val, double prec) const
{
  return ptr_->isInM(val, prec);
}

void ConstManifold::forceOnM(RefVec out, const ConstRefVec& in) const
{
  return ptr_->forceOnM(out, in);
}

Point ConstManifold::getZero() const { return ptr_->getZero(); }

Point ConstManifold::createRandomPoint(double coeff) const
{
  return ptr_->createRandomPoint(coeff);
}

void ConstManifold::createRandomPoint(RefVec out, double coeff) const
{
  return ptr_->createRandomPoint(out, coeff);
}

Index ConstManifold::dim() const { return ptr_->dim(); }

Index ConstManifold::tangentDim() const { return ptr_->tangentDim(); }

Index ConstManifold::representationDim() const
{
  return ptr_->representationDim();
}

void ConstManifold::display(const std::string& prefix) const
{
  return ptr_->display(prefix);
}

void ConstManifold::retractation(RefVec out, const ConstRefVec& x,
                                 const ConstRefVec& v) const
{
  return ptr_->retractation(out, x, v);
}

void ConstManifold::pseudoLog(RefVec out, const ConstRefVec& x,
                              const ConstRefVec& y) const
{
  return ptr_->pseudoLog(out, x, y);
}

void ConstManifold::pseudoLog0(RefVec out, const ConstRefVec& x) const
{
  return ptr_->pseudoLog0(out, x);
}

Eigen::MatrixXd ConstManifold::diffRetractation(const ConstRefVec& x) const
{
  return ptr_->diffRetractation(x);
}

void ConstManifold::applyDiffRetractation(RefMat out, const ConstRefMat& in,
                                          const ConstRefVec& x) const
{
  return ptr_->applyDiffRetractation(out, in, x);
}

Eigen::MatrixXd ConstManifold::diffPseudoLog0(const ConstRefVec& x) const
{
  return ptr_->diffPseudoLog0(x);
}

void ConstManifold::applyDiffPseudoLog0(RefMat out, const ConstRefMat& in,
                                        const ConstRefVec& x) const
{
  return ptr_->applyDiffPseudoLog0(out, in, x);
}

void ConstManifold::applyTransport(RefMat out, const ConstRefMat& in,
                                   const ConstRefVec& x,
                                   const ConstRefVec& v) const
{
  return ptr_->applyTransport(out, in, x, v);
}

void ConstManifold::applyInvTransport(RefMat out, const ConstRefMat& in,
                                      const ConstRefVec& x,
                                      const ConstRefVec& v) const
{
  return ptr_->applyInvTransport(out, in, x, v);
}

void ConstManifold::applyInvTransportOnTheRight(RefMat out,
                                                const ConstRefMat& in,
                                                const ConstRefVec& x,
                                                const ConstRefVec& v) const
{
  return ptr_->applyInvTransportOnTheRight(out, in, x, v);
}

void ConstManifold::setZero(RefVec out) const { return ptr_->setZero(out); }

void ConstManifold::tangentConstraint(RefMat out, const ConstRefVec& x) const
{
  return ptr_->tangentConstraint(out, x);
}

bool ConstManifold::isInTxM(const ConstRefVec& x, const ConstRefVec& v,
                            const double& prec) const
{
  return ptr_->isInTxM(x, v, prec);
}

void ConstManifold::forceOnTxM(RefVec out, const ConstRefVec& in,
                               const ConstRefVec& x) const
{
  return ptr_->forceOnTxM(out, in, x);
}

void ConstManifold::getIdentityOnTxM(RefMat out, const ConstRefVec& x) const
{
  return ptr_->getIdentityOnTxM(out, x);
}

void ConstManifold::limitMap(RefVec out) const { return ptr_->limitMap(out); }

void ConstManifold::getTypicalMagnitude(RefVec out) const
{
  return ptr_->getTypicalMagnitude(out);
}

Eigen::VectorXd ConstManifold::getTypicalMagnitude() const
{
  return ptr_->getTypicalMagnitude();
}

void ConstManifold::getTrustMagnitude(RefVec out) const
{
  return ptr_->getTrustMagnitude(out);
}

Eigen::VectorXd ConstManifold::getTrustMagnitude() const
{
  return ptr_->getTrustMagnitude();
}

void ConstManifold::lock() const { return ptr_->lock(); }

bool ConstManifold::isLocked() const { return ptr_->isLocked(); }

const std::string& ConstManifold::name() const { return ptr_->name(); }

long ConstManifold::getInstanceId() const { return ptr_->getInstanceId(); }

bool ConstManifold::isSameType(const ConstManifold& other) const
{
  return ptr_->isSameType(*(other.ptr_));
}

bool ConstManifold::isElementary() const { return ptr_->isElementary(); }

size_t ConstManifold::numberOfSubmanifolds() const
{
  return ptr_->numberOfSubmanifolds();
}

ConstManifold ConstManifold::operator()(size_t i) const
{
  return ConstManifold(
      std::const_pointer_cast<Manifold_Base>(ptr_->operator()(i)));
}

long ConstManifold::getTypeId() const { return ptr_->getTypeId(); }

std::string ConstManifold::toString(const ConstRefVec& val,
                                    const std::string& prefix,
                                    const Eigen::IOFormat& fmt) const
{
  return ptr_->toString(val, prefix, fmt);
}

/********************************
 *  NON-CONST MANIFOLD METHODS  *
 ********************************/
Manifold::Manifold(std::shared_ptr<Manifold_Base> p) : ConstManifold(p) {}

std::shared_ptr<Manifold_Base> Manifold::getNonConstPtr() { return ptr_; }

void Manifold::setTypicalMagnitude(const ConstRefVec& mag)
{
  return ptr_->setTypicalMagnitude(mag);
}

void Manifold::setTypicalMagnitude(const double& mag)
{
  return ptr_->setTypicalMagnitude(mag);
}

void Manifold::setTrustMagnitude(const ConstRefVec& mag)
{
  return ptr_->setTrustMagnitude(mag);
}

void Manifold::setTrustMagnitude(const double& mag)
{
  return ptr_->setTrustMagnitude(mag);
}

std::string& Manifold::name() { return ptr_->name(); }
}
