#include <iostream>
#include <Eigen/Core>

#include <pgsolver/manifolds/Manifold.h>
#include <pgsolver/manifolds/Point.h>

#include <pgsolver/solver/SolverOptions.h>
#include <pgsolver/solver/LineSearcher.h>

namespace pgs
{
  double LineSearcher::LineSearch(Problem& p, Filter& filter_)
  {
    double alpha = 1.0;
    std::cout << p.x() << std::endl;
    std::cout << filter_.size() << std::endl;
    return alpha;
  }
}
