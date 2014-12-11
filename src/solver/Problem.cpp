#include <pgsolver/solver/Problem.h>

namespace pgs
{
  Problem::Problem(const Manifold& manifold)
    :M_(manifold),
     x_(M_.getIdentity()),
     p_(Eigen::VectorXd::Zero(manifold.dim()))
  {
    std::cout << "Built a Problem" << std::endl;
    std::cout << "x = " << x_ << std::endl;
    std::cout << "p = " << p_ << std::endl;
  }
}
