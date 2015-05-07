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

#ifndef _MANIFOLD_VALID_OBJECT_H_
#define _MANIFOLD_VALID_OBJECT_H_

#include <cstdio>
#include <assert.h>

namespace mnf
{
  class ValidManifold
  {
  private:
    static const int FLAG = 1736274519; //a random number
  public:
    ValidManifold() : flag_(FLAG) {}
    ~ValidManifold() { flag_ = -271828182; }

#ifdef NDEBUG
    bool isValid() const { return true; }

    bool seeMessageAbove() const
    {
      return true;
    }
#else
    bool isValid() const { return flag_ == FLAG; }

    bool seeMessageAbove() const
    {
#  ifndef MNF_ASSERT_THROW
      printf("It appears that you're trying to call a method from a manifold that does not exist \
             anymore. Possible cause: the manifold was statically created in a scope (e.g. a function) which was \
             left since then. For Stan: https://www.youtube.com/watch?v=Yy8MUnlT9Oo\n"); 
#  endif
        return false;
    }
#endif
  private:
    int flag_;
  };
}

#endif //_MANIFOLD_VALID_OBJECT_H_
