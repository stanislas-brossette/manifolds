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
Manifold::Manifold(std::shared_ptr<internal::Manifold_Base> m) : ptr_(m) {std::cout << "Manifold::Manifold(std::shared_ptr<internal::Manifold_Base> m)" << std::endl;}

Manifold::~Manifold() {std::cout << "Manifold::~Manifold()" << std::endl;}

Manifold& Manifold::operator()(size_t)
{std::cout << "Manifold& Manifold::operator()" << std::endl;
  return *this;
}

const Manifold& Manifold::operator()(size_t) const
{std::cout << "const Manifold& Manifold::operator()" << std::endl;
  return *this;
}

Manifold Manifold::shallowCopy()
{std::cout << "Manifold Manifold::shallowCopy()" << std::endl;
  return (Manifold(ptr_));
}

Manifold Manifold::deepCopy() const
{std::cout << "Manifold Manifold::deepCopy()" << std::endl;
  return (Manifold(ptr_->clone()));
}

std::shared_ptr<const internal::Manifold_Base> Manifold::ptr() const
{std::cout << "std::shared_ptr<const internal::Manifold_Base> Manifold::ptr()" << std::endl;
  return std::const_pointer_cast<const internal::Manifold_Base>(ptr_);
}

Point Manifold::createPoint() const 
{std::cout << "Point Manifold::createPoint()" << std::endl; 
  lock();
  return Point(*this);
}

Point Manifold::createPoint(const ConstRefVec& val) const
{std::cout << "Point Manifold::createPoint(const ConstRefVec& val)" << std::endl;
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

bool Manifold::isInM(const Eigen::VectorXd& val, double prec) const
{
  return ptr_->isInM(val, prec);
}

void Manifold::forceOnM(RefVec out, const ConstRefVec& in) const
{
  return ptr_->forceOnM(out, in);
}

Point Manifold::getZero() const 
{ 
  Point p = createPoint();
  setZero(p.value());
  return p; 
}

Point Manifold::createRandomPoint(double coeff) const
{
  Point p = createPoint();
  createRandomPoint(p.value(), coeff);
  return p; 
}

void Manifold::createRandomPoint(RefVec out, double coeff) const
{
  return ptr_->createRandomPoint(out, coeff);
}

Index Manifold::dim() const { return ptr_->dim(); }

Index Manifold::tangentDim() const { return ptr_->tangentDim(); }

Index Manifold::representationDim() const
{
  return ptr_->representationDim();
}

void Manifold::display(const std::string& prefix) const
{
  return ptr_->display(prefix);
}

void Manifold::retractation(RefVec out, const ConstRefVec& x,
                                 const ConstRefVec& v) const
{
  return ptr_->retractation(out, x, v);
}

void Manifold::pseudoLog(RefVec out, const ConstRefVec& x,
                              const ConstRefVec& y) const
{
  return ptr_->pseudoLog(out, x, y);
}

void Manifold::pseudoLog0(RefVec out, const ConstRefVec& x) const
{
  return ptr_->pseudoLog0(out, x);
}

Eigen::MatrixXd Manifold::diffRetractation(const ConstRefVec& x) const
{
  return ptr_->diffRetractation(x);
}

void Manifold::applyDiffRetractation(RefMat out, const ConstRefMat& in,
                                          const ConstRefVec& x) const
{
  return ptr_->applyDiffRetractation(out, in, x);
}

Eigen::MatrixXd Manifold::diffPseudoLog0(const ConstRefVec& x) const
{
  return ptr_->diffPseudoLog0(x);
}

void Manifold::applyDiffPseudoLog0(RefMat out, const ConstRefMat& in,
                                        const ConstRefVec& x) const
{
  return ptr_->applyDiffPseudoLog0(out, in, x);
}

void Manifold::applyTransport(RefMat out, const ConstRefMat& in,
                                   const ConstRefVec& x,
                                   const ConstRefVec& v) const
{
  return ptr_->applyTransport(out, in, x, v);
}

void Manifold::applyInvTransport(RefMat out, const ConstRefMat& in,
                                      const ConstRefVec& x,
                                      const ConstRefVec& v) const
{
  return ptr_->applyInvTransport(out, in, x, v);
}

void Manifold::applyInvTransportOnTheRight(RefMat out,
                                                const ConstRefMat& in,
                                                const ConstRefVec& x,
                                                const ConstRefVec& v) const
{
  return ptr_->applyInvTransportOnTheRight(out, in, x, v);
}

void Manifold::setZero(RefVec out) const { return ptr_->setZero(out); }

void Manifold::tangentConstraint(RefMat out, const ConstRefVec& x) const
{
  return ptr_->tangentConstraint(out, x);
}

bool Manifold::isInTxM(const ConstRefVec& x, const ConstRefVec& v,
                            const double& prec) const
{
  return ptr_->isInTxM(x, v, prec);
}

void Manifold::forceOnTxM(RefVec out, const ConstRefVec& in,
                               const ConstRefVec& x) const
{
  return ptr_->forceOnTxM(out, in, x);
}

void Manifold::getIdentityOnTxM(RefMat out, const ConstRefVec& x) const
{
  return ptr_->getIdentityOnTxM(out, x);
}

void Manifold::limitMap(RefVec out) const { return ptr_->limitMap(out); }

void Manifold::getTypicalMagnitude(RefVec out) const
{
  return ptr_->getTypicalMagnitude(out);
}

Eigen::VectorXd Manifold::getTypicalMagnitude() const
{
  return ptr_->getTypicalMagnitude();
}

void Manifold::getTrustMagnitude(RefVec out) const
{
  return ptr_->getTrustMagnitude(out);
}

Eigen::VectorXd Manifold::getTrustMagnitude() const
{
  return ptr_->getTrustMagnitude();
}

void Manifold::lock() const { return ptr_->lock(); }

bool Manifold::isLocked() const { return ptr_->isLocked(); }

const std::string& Manifold::name() const { return ptr_->name(); }
std::string& Manifold::name() { return ptr_->name(); }

long Manifold::getInstanceId() const { return ptr_->getInstanceId(); }

bool Manifold::isSameType(const Manifold& other) const
{
  return ptr_->isSameType(*(other.ptr_));
}

bool Manifold::isElementary() const { return ptr_->isElementary(); }

size_t Manifold::numberOfSubmanifolds() const
{
  return ptr_->numberOfSubmanifolds();
}

long Manifold::getTypeId() const { return ptr_->getTypeId(); }

std::string Manifold::toString(const ConstRefVec& val,
                                    const std::string& prefix,
                                    const Eigen::IOFormat& fmt) const
{
  return ptr_->toString(val, prefix, fmt);
}

std::shared_ptr<internal::Manifold_Base> Manifold::getNonConstPtr() {std::cout << "std::shared_ptr<internal::Manifold_Base> Manifold::getNonConstPtr()" << std::endl; return ptr_; }

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

}
