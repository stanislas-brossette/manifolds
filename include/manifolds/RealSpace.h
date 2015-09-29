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

#include <iostream>
#include <manifolds/defs.h>
#include <manifolds/Manifold.h>
#include <manifolds/RealSpace_Base.h>
#include <manifolds/utils.h>

namespace mnf
{
  /// \brief Manifold representing the space of real numbers of dimension n
  /// \f$\mathbb{R}^n\f$
  class MANIFOLDS_API RealSpace: public Manifold
  {
  public:
    /// \brief Constructor
    /// \param n the dimension of the realspace \f$\mathbb{R}^n\f$
    RealSpace(Index n);
    RealSpace(Index n, double magnitude);
    RealSpace(Index n, const ConstRefVec& magnitude);

    virtual size_t numberOfSubmanifolds() const;
    virtual Manifold operator()(size_t i) const;

    virtual std::string toString(const ConstRefVec& val, const std::string& prefix = "", const Eigen::IOFormat& fmt = mnf::defaultFormat) const;
    virtual void getTypicalMagnitude_(RefVec out) const;
    void setTypicalMagnitude(double magnitude);
    void setTypicalMagnitude(const ConstRefVec& out);
    virtual void getTrustMagnitude_(RefVec out) const;
    void setTrustMagnitude(const double& magnitude);
    void setTrustMagnitude(const ConstRefVec& out);
    virtual bool isElementary() const;
    virtual long getTypeId() const;

  };
}

