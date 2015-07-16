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

#ifndef _WIN32
#define BOOST_TEST_MODULE Manifolds
#endif

#include <boost/test/unit_test.hpp>

#include <manifolds/defs.h>
#include <manifolds/utils.h>
#include <manifolds/Point.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>

using namespace mnf;

// Asserts with exception throwing are performed only in debug
#ifdef NDEBUG
# define CHECK_THROW_IN_DEBUG(expression, exception) ((void)0)
#else
# define CHECK_THROW_IN_DEBUG(expression, exception) BOOST_CHECK_THROW(expression, exception)
#endif
#ifdef NDEBUG
# define CHECK_NO_THROW_IN_DEBUG(expression, exception) ((void)0)
#else
# define CHECK_NO_THROW_IN_DEBUG(expression, exception) BOOST_CHECK_NO_THROW(expression)
#endif


// This function is not bugged anymore; since the manifolds stored
// in the CartesianProduct are copies of the original manifold, it
// will create the Point without any issues.
CartesianProduct* buildProduct()
{
  RealSpace R3(3);
  RealSpace R2(2);
  CartesianProduct* M = new CartesianProduct(R3, R2);
  return M;
}

BOOST_AUTO_TEST_CASE(ManifoldIsValid)
{
  CHECK_NO_THROW_IN_DEBUG(buildProduct()->getZero(), mnf::mnf_exception);
}

Point createR3Point()
{
  RealSpace R3(3);
  return R3.createPoint();
}

BOOST_AUTO_TEST_CASE(ManifoldRefCounter)
{
  CHECK_THROW_IN_DEBUG(createR3Point(), mnf::mnf_exception);
}
