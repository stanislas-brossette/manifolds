#include <pgsolver/solver/Solver.h>
#include <pgsolver/solver/ConstraintManager.h>

namespace pgs
{
  Solver::Solver()
  {
    std::cout << "New Solver" << std::endl;
  }

  Results Solver::solve(Problem& problem, Point& x0)
  {
    cstrMngr_.init(problem);
    initSolver(problem, x0);
    problem.printState();
    //int iter = 0;
    //int maxIter = 1000;
    //double epsilon_P = 1e-6;
    //double epsilon_D = 1e-6;

    //Setup the matrix to hold the 

    std::cout << opt_.epsilon_P << std::endl;

    return Results({ x0, CONVERGE, {} });
  }

  void Solver::initSolver(Problem& problem, Point& x0)
  {
    problem.setX(x0);
    opt_.maxIter = 10000;
    opt_.epsilon_P = 1e-6;
    opt_.epsilon_D = 1e-2;
  }

}
