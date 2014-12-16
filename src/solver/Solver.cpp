#include <pgsolver/solver/Solver.h>

namespace pgs
{
  void Solver::solve(Problem& problem, Point& x0)
  {
    problem.setX(x0);
    problem.printState();
  }

}
