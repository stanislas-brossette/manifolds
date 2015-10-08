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

#ifndef _MANIFOLDS_REAL_SPACE_H_
#define _MANIFOLDS_REAL_SPACE_H_

#include <iostream>
#include <manifolds/defs.h>
#include <manifolds/Manifold.h>
#include <manifolds/utils.h>

namespace mnf
{
/// \brief Manifold representing the space of real numbers of dimension n
/// \f$\mathbb{R}^n\f$
class MANIFOLDS_API RealSpace : public Manifold
{
 public:
  /// \brief Constructor
  /// \param n the dimension of the realspace \f$\mathbb{R}^n\f$
  RealSpace(Index n);
  RealSpace(Index n, double magnitude);
  RealSpace(Index n, const ConstRefVec& magnitude);

  virtual size_t numberOfSubmanifolds() const;
  virtual const Manifold& operator()(size_t i) const;

  virtual std::string toString(
      const ConstRefVec& val, const std::string& prefix = "",
      const Eigen::IOFormat& fmt = mnf::defaultFormat) const;
  virtual void getTypicalMagnitude_(RefVec out) const;
  void setTypicalMagnitude(double magnitude);
  void setTypicalMagnitude(const ConstRefVec& out);
  virtual void getTrustMagnitude_(RefVec out) const;
  void setTrustMagnitude(const double& magnitude);
  void setTrustMagnitude(const ConstRefVec& out);
  virtual bool isElementary() const;
  virtual long getTypeId() const;

 protected:
  // map operations
  virtual bool isInM_(const Eigen::VectorXd&, double prec) const;
  virtual void forceOnM_(RefVec out, const ConstRefVec& in) const;
  virtual void getIdentityOnTxM_(RefMat out, const ConstRefVec& x) const;
  virtual void createRandomPoint_(RefVec out, double coeff) const;
  virtual void retractation_(RefVec out, const ConstRefVec& x,
                             const ConstRefVec& v) const;
  virtual void pseudoLog_(RefVec out, const ConstRefVec& x,
                          const ConstRefVec& y) const;
  virtual void pseudoLog0_(RefVec out, const ConstRefVec& x) const;
  virtual void setZero_(RefVec out) const;
  virtual Eigen::MatrixXd diffRetractation_(const ConstRefVec& x) const;
  virtual void applyDiffRetractation_(RefMat out, const ConstRefMat& in,
                                      const ConstRefVec& x) const;
  virtual Eigen::MatrixXd diffPseudoLog0_(const ConstRefVec& x) const;
  virtual void applyDiffPseudoLog0_(RefMat out, const ConstRefMat& in,
                                    const ConstRefVec& x) const;
  virtual void applyTransport_(RefMat out, const ConstRefMat& in,
                               const ConstRefVec& x,
                               const ConstRefVec& v) const;
  virtual void applyInvTransport_(RefMat out, const ConstRefMat& in,
                                  const ConstRefVec& x,
                                  const ConstRefVec& v) const;
  virtual void applyInvTransportOnTheRight_(RefMat out, const ConstRefMat& in,
                                            const ConstRefVec& x,
                                            const ConstRefVec& v) const;

  virtual void tangentConstraint_(RefMat out, const ConstRefVec& x) const;
  virtual bool isInTxM_(const ConstRefVec& x, const ConstRefVec& v,
                        const double& prec) const;
  virtual void forceOnTxM_(RefVec out, const ConstRefVec& in,
                           const ConstRefVec& x) const;
  virtual void limitMap_(RefVec out) const;

  virtual Manifold_ptr getNewCopy() const;

 private:
  Eigen::VectorXd typicalMagnitude_;
  Eigen::VectorXd trustMagnitude_;
};
}

#endif  //_MANIFOLDS_REAL_SPACE_H_
