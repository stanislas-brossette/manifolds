#ifndef _WIN32
#define BOOST_TEST_MODULE PGSolver 
#endif

#include <boost/test/unit_test.hpp>

#include <manifolds/defs.h>
#include <manifolds/utils.h>
#include <manifolds/Point.h>
#include <manifolds/RealSpace.h>
#include <manifolds/CartesianProduct.h>

using namespace pgs;

// Asserts with exception throwing are performed only in debug
#ifdef NDEBUG
# define CHECK_THROW_IN_DEBUG(expression, exception) ((void)0)
#else
# define CHECK_THROW_IN_DEBUG(expression, exception) BOOST_CHECK_THROW(expression, exception)
#endif

// A bugged helper function with some statically created variables that will
// be used outside of the function scope.
CartesianProduct* buildProduct()
{
  RealSpace R3(3);
  RealSpace R2(2);
  CartesianProduct* M = new CartesianProduct(R3, R2);
  return M;
}

BOOST_AUTO_TEST_CASE(ManifoldIsValid)
{
  CartesianProduct* M = buildProduct();
  CHECK_THROW_IN_DEBUG(M->getIdentity(), pgs::pgs_exception);
}

Point createR3Point()
{
  RealSpace R3(3);
  return R3.createPoint();
}

BOOST_AUTO_TEST_CASE(ManifoldRefCounter)
{
  CHECK_THROW_IN_DEBUG(createR3Point(), pgs::pgs_exception);
}
