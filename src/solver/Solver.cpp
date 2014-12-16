#include <pgsolver/solver/Solver.h>

namespace pgs
{
  Solver::Solver()
  {
  }

  Results Solver::solve(Problem& problem, Point& x0)
  {
    problem.setX(x0);
    problem.printState();

    return Results({ x0, CONVERGE, {} });
  }

}
