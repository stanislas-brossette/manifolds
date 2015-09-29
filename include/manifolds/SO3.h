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
#include <manifolds/Manifold.h>
#include <manifolds/ReusableTemporaryMap.h>
#include <manifolds/utils.h>

namespace mnf
{
  template<typename Map>
  class SO3_Base<Map>;
  /// \brief Manifold representing the space of 3-dimensional rotations, also
  /// known as SO(3). It is templated by its map
  template<typename Map>
  class SO3: public Manifold
  {
  public:
    SO3();
    SO3(double magnitude);
    SO3(const ConstRefVec& magnitude);
    virtual size_t numberOfSubmanifolds() const;
    virtual const Manifold& operator()(size_t i) const;
    virtual std::string toString(const ConstRefVec& val, const std::string& prefix = "", const Eigen::IOFormat& fmt = mnf::defaultFormat) const;
    virtual void getTypicalMagnitude_(RefVec out) const;
    void setTypicalMagnitude(double magnitude);
    void setTypicalMagnitude(const ConstRefVec& out);
    virtual void getTrustMagnitude_(RefVec out) const;
    void setTrustMagnitude(const double& magnitude);
    void setTrustMagnitude(const ConstRefVec& out);
    virtual void createRandomPoint_(RefVec out, double coeff) const;
    virtual bool isElementary() const;
    virtual long getTypeId() const;
  };

  //Implementations of the methods
  template<typename Map>
  inline SO3<Map>::SO3()
    : Manifold(std::make_shared<SO3_Base<Map> >())
  {
  }
  template<typename Map>
  inline SO3<Map>::SO3(double magnitude)
    : Manifold(std::make_shared<SO3_Base<Map> >(magnitude))
  {
  }
  template<typename Map>
  inline SO3<Map>::SO3(const ConstRefVec& magnitude)
    : Manifold(std::make_shared<SO3_Base<Map> >(magnitude))
  {
  }

  template<typename Map>
  inline size_t SO3<Map>::numberOfSubmanifolds() const
  {
    return std::static_pointer_cast<SO3_Base<Map> >(manifoldBase_)->numberOfSubmanifolds();
  }

  template<typename Map>
  inline bool SO3<Map>::isElementary() const
  {
    return std::static_pointer_cast<SO3_Base<Map> >(manifoldBase_)->isElementary();
  }

  template<typename Map>
  inline const Manifold& SO3<Map>::operator()(size_t i) const
  {
    return std::static_pointer_cast<SO3_Base<Map> >(manifoldBase_)->operator()(i);
  }

  template<typename Map>
  inline std::string SO3<Map>::toString(const ConstRefVec& val, const std::string& prefix, const Eigen::IOFormat& fmt) const
  {
    return std::static_pointer_cast<SO3_Base<Map> >(manifoldBase_)->toString(val, prefix, fmt);
  }

  template<typename Map>
  void SO3<Map>::createRandomPoint_(RefVec out, double coeff) const
  {
    return std::static_pointer_cast<SO3_Base<Map> >(manifoldBase_)->createRandomPoint_(out, coeff);
  }

  template<typename Map>
  void SO3<Map>::getTypicalMagnitude_(RefVec out) const
  {
    return std::static_pointer_cast<SO3_Base<Map> >(manifoldBase_)->getTypicalMagnitude_(out);
  }

  template<typename Map>
  void SO3<Map>::setTypicalMagnitude(double magnitude)
  {
    return std::static_pointer_cast<SO3_Base<Map> >(manifoldBase_)->setTypicalMagnitude(magnitude);
  }

  template<typename Map>
  void SO3<Map>::setTypicalMagnitude(const ConstRefVec& out)
  {
    return std::static_pointer_cast<SO3_Base<Map> >(manifoldBase_)->setTypicalMagnitude(out);
  }

  template<typename Map>
  void SO3<Map>::getTrustMagnitude_(RefVec out) const
  {
    return std::static_pointer_cast<SO3_Base<Map> >(manifoldBase_)->getTrustMagnitude_(out);
  }

  template<typename Map>
  void SO3<Map>::setTrustMagnitude(const double& magnitude)
  {
    return std::static_pointer_cast<SO3_Base<Map> >(manifoldBase_)->setTrustMagnitude(magnitude);
  }

  template<typename Map>
  void SO3<Map>::setTrustMagnitude(const ConstRefVec& out)
  {
    return std::static_pointer_cast<SO3_Base<Map> >(manifoldBase_)->setTrustMagnitude(out);
  }

  template<typename Map>
  long SO3<Map>::getTypeId() const
  {
    return std::static_pointer_cast<SO3_Base<Map> >(manifoldBase_)->getTypeId();
  }
}
