#include <pgsolver/solver/Solver.cpp>

namespace pgs
{
  Solver::Solver()
  {
  }

  void Solver::solve(Problem& problem, Point& x0)
  {
    problem.setX(x0);
    problem.printState();
  }

}
