#ifndef _MANIFOLDS_CARTESIAN_POWER_H_
#define _MANIFOLDS_CARTESIAN_POWER_H_

#include <manifolds/defs.h>
#include <manifolds/CartesianProduct.h>

namespace pgs
{
  /// \brief Manifold representing the cartesian product of n times the same manifold
  class MANIFOLDS_API CartesianPower : public CartesianProduct
  {
  public:
    /// \brief Constructor of the \f$ M^n m2\f$
    CartesianPower(const Manifold& M, int n);
  };

}


#endif //_MANIFOLDS_CARTESIAN_POWER_H_