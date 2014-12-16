#include <pgsolver/solver/Solver.h>

namespace pgs
{
  Solver::Solver()
  {
    std::cout << "New Solver" << std::endl;
  }

  Results Solver::solve(Problem& problem, Point& x0)
  {
    problem.setX(x0);
    problem.printState();

    Results r;
    return r;
  }

}
