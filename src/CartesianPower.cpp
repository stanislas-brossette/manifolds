#include <manifolds/CartesianPower.h>

namespace pgs
{
  CartesianPower::CartesianPower(const Manifold& M, int n)
    : CartesianProduct()
  {
    for (int i = 0; i < n; ++i)
      multiply(M);
  }
}


