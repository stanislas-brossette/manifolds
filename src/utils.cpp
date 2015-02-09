#include <manifolds/defs.h>
#include <manifolds/utils.h>

#include <Eigen/Core>

namespace utils
{
  bool set_is_malloc_allowed(bool allow)
  {
#ifdef EIGEN_RUNTIME_NO_MALLOC
    return Eigen::internal::set_is_malloc_allowed(allow);
#else
    eigen_assert(false && "you can't call this function if EigenQP was compiled without the flag EIGEN_RUNTIME_NO_MALLOC");
    return true;
#endif
  }
}
