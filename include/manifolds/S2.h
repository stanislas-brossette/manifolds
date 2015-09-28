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

#ifndef _MANIFOLDS_S2_H_
#define _MANIFOLDS_S2_H_
#include <stdexcept>

#include <manifolds/defs.h>
#include <manifolds/Manifold.h>
#include <manifolds/utils.h>

namespace mnf
{
  class S2_Base;

  /// \brief Manifold representing the 3-dimensional Sphere, also
  /// known as S2.
  /// All the equations in this class are provided by Manopt
  class MANIFOLDS_API S2 : public Manifold
  {
  public:
    S2();
    S2(double trustMagnitude);
    S2(const ConstRefVec& trustMagnitude);

    virtual size_t numberOfSubmanifolds() const;
    virtual const Manifold& operator()(size_t i) const;

    virtual void createRandomPoint_(RefVec out, double coeff) const;

    virtual std::string toString(const ConstRefVec& val, const std::string& prefix = "", const Eigen::IOFormat& fmt = mnf::defaultFormat) const;
    virtual void getTypicalMagnitude_(RefVec out) const;
    void setTypicalMagnitude(double magnitude);
    void setTypicalMagnitude(const ConstRefVec& out);
    virtual void getTrustMagnitude_(RefVec out) const;
    void setTrustMagnitude(const double& magnitude);
    void setTrustMagnitude(const ConstRefVec& out);
    void logarithm (RefVec out, const ConstRefVec& x, const ConstRefVec& y) const;
    double distance (const ConstRefVec& x, const ConstRefVec& y) const;
    /// \brief projects each row of \ in on TxM
    void projRows (RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    /// \brief projects each cols of \ in on TxM
    void projCols (RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    /// \brief projects vector \ in on TxM
    void projVec (RefVec out, const ConstRefVec& in, const ConstRefVec& x) const;
    void rand(RefVec out) const;
    void randVec(RefVec out, const ConstRefVec& x) const;
    Eigen::Vector3d randVec(const ConstRefVec& x) const;
    virtual bool isElementary() const;
    virtual long getTypeId() const;

  };
}
#endif //_MANIFOLDS_S2_H_
