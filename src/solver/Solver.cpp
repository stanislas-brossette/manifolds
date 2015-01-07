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
    problem.setX(x0);
    cstrMngr_.init(problem);
    initSolver(problem);
    std::cout << "Problem with Linear cstr of Dim: " << cstrMngr_.totalDimLin()<< std::endl;
    std::cout << "And NonLinear cstr of Dim: " << cstrMngr_.totalDimNonLin()<< std::endl;
    updateAllProblemData(problem);
    printStatus();

    Eigen::Vector3d z1(1,1,1);
    problem.setX(x0+z1);
    updateAllProblemData(problem);
    printStatus();
    //int iter = 0;
    //int maxIter = 1000;
    //double epsilon_P = 1e-6;
    //double epsilon_D = 1e-6;

    return Results({ x0, CONVERGE, {} });
  }

  void Solver::printStatus()
  {
    std::cout << "================================================================="<< std::endl;
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    std::cout << "current x = " << problem_->x() << std::endl;
    std::cout << "current z = " << z_.transpose().format(CleanFmt) << std::endl;
    probEval_.print();
    std::cout << "================================================================="<< std::endl;
  }

  void Solver::initSolver(Problem& problem)
  {
    problem_ = &problem;
    opt_.maxIter = 10000;
    opt_.epsilon_P = 1e-6;
    opt_.epsilon_D = 1e-2;
    probEval_.diffObj.resize(1, problem.M().dim());
    probEval_.tangentLB.resize(problem.M().dim());
    probEval_.tangentUB.resize(problem.M().dim());

    probEval_.linCstr.resize(cstrMngr_.totalDimLin());
    probEval_.diffLinCstr.resize(cstrMngr_.totalDimLin(), problem.M().dim());
    probEval_.linCstrLB.resize(cstrMngr_.totalDimLin());
    probEval_.linCstrUB.resize(cstrMngr_.totalDimLin());

    probEval_.nonLinCstr.resize(cstrMngr_.totalDimNonLin());
    probEval_.diffNonLinCstr.resize(cstrMngr_.totalDimNonLin(), problem.M().dim());
    probEval_.nonLinCstrLB.resize(cstrMngr_.totalDimNonLin());
    probEval_.nonLinCstrUB.resize(cstrMngr_.totalDimNonLin());

    probEval_.Hessian.resize(problem.M().dim(), problem.M().dim());

    z_.resize(problem.M().dim());
    lagMultLin_.resize(cstrMngr_.totalDimLin());
    lagMultNonLin_.resize(cstrMngr_.totalDimNonLin());
  }

  void Solver::updateAllProblemData(Problem& p)
  {
    p.evalObj(probEval_.obj); 
    p.evalObjDiff(probEval_.diffObj); 
    p.getTangentLB(probEval_.tangentLB);
    p.getTangentUB(probEval_.tangentUB);

    for (size_t i = 0; i<p.numberOfCstr(); ++i)
    {
      // for each constraint we fill the correct lines of the matrices using
      // getView
      p.evalLinCstr(cstrMngr_.getViewLin(probEval_.linCstr,i),i);
      p.evalLinCstrDiff(cstrMngr_.getViewLin(probEval_.diffLinCstr,i),i);
      p.getLinCstrLB(cstrMngr_.getViewLin(probEval_.linCstrLB,i),i);

      p.evalNonLinCstr(cstrMngr_.getViewNonLin(probEval_.nonLinCstr,i),i);
      p.evalNonLinCstrDiff(cstrMngr_.getViewNonLin(probEval_.diffNonLinCstr,i),i);
    }
  }
}
