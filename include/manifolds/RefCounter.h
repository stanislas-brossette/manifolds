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

#include <manifolds/defs.h>
#include <manifolds/Point.h>
#include <manifolds/mnf_assert.h>

// The operator version of noexcept is available starting Visual Studio 2015
// https://msdn.microsoft.com/en-us/library/dn956976.aspx
#if defined(_MSC_FULL_VER) && _MSC_VER < 1900
#define NOEXCEPT(x)
#else
#define NOEXCEPT(x) noexcept(x)
#endif

namespace mnf
{
  /// \brief object containing a counter and that cannot be destroyed if the
  /// counter is not at 0.
  class RefCounter
  {
    public:
      RefCounter()
        :count_(0)
      {
      }

  RefCounter(const RefCounter&)
    :count_(0)
      {
      }

      // In C++ 11, by default, a destructor cannot launch an exception, setting noexcet to false allows it.
      ~RefCounter() NOEXCEPT(false)
      {
        mnf_assert(count_ == 0 && "You cannot destroy this manifold because some points still depend on it");
      }

    protected:
      void incrementRefCounter() const
      {
        count_++;
      }
      void decrementRefCounter() const
      {
        mnf_assert(count_>0 && "You cannot decrement when no point exist");
        count_--;
      }

    private:
      mutable int count_;
      friend void ConstSubPoint::registerPoint();
      friend void ConstSubPoint::unregisterPoint();
  };
}

