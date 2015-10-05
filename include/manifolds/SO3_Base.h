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

#pragma once
#define _USE_MATH_DEFINES
#include <math.h>

#include <manifolds/defs.h>
#include <manifolds/Manifold_Base.h>
#include <manifolds/SO3.h>
#include <manifolds/ReusableTemporaryMap.h>
#include <manifolds/utils.h>

namespace mnf
{
/// \brief Manifold representing the space of 3-dimensional rotations, also
/// known as SO(3). It is templated by its map
template <typename Map>
class MANIFOLDS_API SO3_Base : public Manifold_Base
{
  template <typename>
  friend class SO3;

 private:
  SO3_Base();
  SO3_Base(double magnitude);
  SO3_Base(const ConstRefVec& magnitude);

  virtual SO3_Base copy() const;
  virtual std::shared_ptr<Manifold_Base> clone() const;

  virtual size_t numberOfSubmanifolds() const;
  virtual std::shared_ptr<const Manifold_Base> operator()(size_t i) const;
  virtual std::string toString(
      const ConstRefVec& val, const std::string& prefix = "",
      const Eigen::IOFormat& fmt = mnf::defaultFormat) const;
  virtual void createRandomPoint_(RefVec out, double coeff) const;
  virtual bool isElementary() const;
  virtual long getTypeId() const;

 protected:
  // map operations
  virtual bool isInM_(const Eigen::VectorXd& val, double prec) const;
  virtual void forceOnM_(RefVec out, const ConstRefVec& in) const;
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
  virtual void getIdentityOnTxM_(RefMat out, const ConstRefVec& x) const;
  virtual void limitMap_(RefVec out) const;
  mutable ReusableTemporaryMap bufferMap_;

 private:
  Eigen::Vector3d typicalMagnitude_;
  Eigen::Vector3d trustMagnitude_;
};

// Implementations of the methods
template <typename Map>
inline SO3_Base<Map>::SO3_Base()
    : Manifold_Base(3, Map::InputDim_, Map::OutputDim_)
{
  name() = "SO3";
  setTypicalMagnitude(Eigen::Vector3d::Constant(M_PI));
  setTrustMagnitude(Eigen::Vector3d::Constant(M_PI));
}
template <typename Map>
inline SO3_Base<Map>::SO3_Base(double magnitude)
    : Manifold_Base(3, Map::InputDim_, Map::OutputDim_)
{
  name() = "SO3";
  setTypicalMagnitude(Eigen::Vector3d::Constant(magnitude));
  setTrustMagnitude(Eigen::Vector3d::Constant(magnitude));
}
template <typename Map>
inline SO3_Base<Map>::SO3_Base(const ConstRefVec& magnitude)
    : Manifold_Base(3, Map::InputDim_, Map::OutputDim_)
{
  mnf_assert(magnitude.size() == 3 && "magnitude on SO3 must be of size 3");
  name() = "SO3";
  setTypicalMagnitude(magnitude);
  setTrustMagnitude(magnitude);
}

template <typename Map>
inline bool SO3_Base<Map>::isInM_(const Eigen::VectorXd& val, double prec) const
{
  return Map::isInM_(val, prec);
}

template <typename Map>
inline void SO3_Base<Map>::forceOnM_(RefVec out, const ConstRefVec& in) const
{
  return Map::forceOnM_(out, in);
}

template <typename Map>
inline SO3_Base<Map> SO3_Base<Map>::copy() const
{
  return SO3_Base<Map>(*this);
}

template <typename Map>
inline std::shared_ptr<Manifold_Base> SO3_Base<Map>::clone() const
{
  return std::make_shared<SO3_Base<Map>>(*this);
}

template <typename Map>
inline size_t SO3_Base<Map>::numberOfSubmanifolds() const
{
  return 1;
}

template <typename Map>
inline bool SO3_Base<Map>::isElementary() const
{
  return true;
}

template <typename Map>
inline std::shared_ptr<const Manifold_Base> SO3_Base<Map>::operator()(
    size_t) const
{
  return shared_from_this();
}

template <typename Map>
inline std::string SO3_Base<Map>::toString(const ConstRefVec& val,
                                           const std::string& prefix,
                                           const Eigen::IOFormat& fmt) const
{
  std::string matPrefix = prefix + '[';
  std::stringstream ss;
  ss << matPrefix
     << (Eigen::Map<const typename Map::DisplayType>(val.data())).format(fmt);
  return ss.str();
}

template <typename Map>
void SO3_Base<Map>::createRandomPoint_(RefVec out, double coeff) const
{
  setZero_(out);
  Eigen::Vector3d v(coeff * Eigen::Vector3d::Random());
  retractation(out, out, v);
}

template <typename Map>
inline void SO3_Base<Map>::retractation_(RefVec out, const ConstRefVec& x,
                                         const ConstRefVec& v) const
{
  Map::retractation_(out, x, v);
}

template <typename Map>
inline void SO3_Base<Map>::pseudoLog_(RefVec out, const ConstRefVec& x,
                                      const ConstRefVec& y) const
{
  Map::pseudoLog_(out, x, y);
}

template <typename Map>
inline void SO3_Base<Map>::pseudoLog0_(RefVec out, const ConstRefVec& x) const
{
  Map::pseudoLog0_(out, x);
}

template <typename Map>
inline void SO3_Base<Map>::setZero_(RefVec out) const
{
  Map::setZero_(out);
}

template <typename Map>
inline Eigen::MatrixXd SO3_Base<Map>::diffRetractation_(
    const ConstRefVec& x) const
{
  return Map::diffRetractation_(x);
}

template <typename Map>
inline void SO3_Base<Map>::applyDiffRetractation_(RefMat out,
                                                  const ConstRefMat& in,
                                                  const ConstRefVec& x) const
{
  Map::applyDiffRetractation_(out, in, x, bufferMap_);
}

template <typename Map>
inline Eigen::MatrixXd SO3_Base<Map>::diffPseudoLog0_(
    const ConstRefVec& x) const
{
  return Map::diffPseudoLog0_(x);
}

template <typename Map>
inline void SO3_Base<Map>::applyDiffPseudoLog0_(RefMat out,
                                                const ConstRefMat& in,
                                                const ConstRefVec& x) const
{
  Map::applyDiffPseudoLog0_(out, in, x, bufferMap_);
}

template <typename Map>
inline void SO3_Base<Map>::applyTransport_(RefMat out, const ConstRefMat& in,
                                           const ConstRefVec& x,
                                           const ConstRefVec& v) const
{
  Map::applyTransport_(out, in, x, v);
}

template <typename Map>
inline void SO3_Base<Map>::applyInvTransport_(RefMat out, const ConstRefMat& in,
                                              const ConstRefVec& x,
                                              const ConstRefVec& v) const
{
  Map::applyInvTransport_(out, in, x, v);
}

template <typename Map>
inline void SO3_Base<Map>::applyInvTransportOnTheRight_(
    RefMat out, const ConstRefMat& in, const ConstRefVec& x,
    const ConstRefVec& v) const
{
  Map::applyInvTransportOnTheRight_(out, in, x, v);
}

template <typename Map>
void SO3_Base<Map>::tangentConstraint_(RefMat out, const ConstRefVec& x) const
{
  Map::tangentConstraint_(out, x);
}

template <typename Map>
bool SO3_Base<Map>::isInTxM_(const ConstRefVec& x, const ConstRefVec& v,
                             const double& prec) const
{
  return Map::isInTxM_(x, v, prec);
}

template <typename Map>
void SO3_Base<Map>::forceOnTxM_(RefVec out, const ConstRefVec& in,
                                const ConstRefVec& x) const
{
  Map::forceOnTxM_(out, in, x);
}

template <typename Map>
void SO3_Base<Map>::getIdentityOnTxM_(RefMat out, const ConstRefVec& x) const
{
  Map::getIdentityOnTxM_(out, x);
}

template <typename Map>
void SO3_Base<Map>::limitMap_(RefVec out) const
{
  double limitExpMap = M_PI / sqrt(3);
  out.setConstant(limitExpMap);
}

template <typename Map>
long SO3_Base<Map>::getTypeId() const
{
  constexpr long typeId = utils::hash::computeHash("SO3", Map::hashName);
  return typeId;
}
}