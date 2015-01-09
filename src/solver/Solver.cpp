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

    z_.setZero();
    lagMult_.linear.setOnes();
    lagMult_.nonLinear.setOnes();
    lagMult_.update(lagMult_.linear, lagMult_.nonLinear);
    lagMult_.print();
    updateAllProblemData(problem);
    printStatus();

    //int iter = 0;
    //int maxIter = 1000;
    //double epsilon_P = 1e-6;
    //double epsilon_D = 1e-6;
    bool converged = convergence(opt_.epsilon_P , opt_.epsilon_P, 
        x0, lagMult_.all, probEval_.allCstr, probEval_.diffLag);
    std::cout << converged << std::endl;
     

    return Results({ x0, CONVERGE, {} });
  }

  void Solver::printStatus()
  {
    std::cout << "================================================================="<< std::endl;
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    std::cout << "current x = " << problem_->x() << std::endl;
    std::cout << "current z = " << z_.transpose().format(CleanFmt) << std::endl;
    std::cout << "current Lagrange mult for Lin Cstr: " << lagMult_.linear.transpose().format(CleanFmt) << std::endl;
    std::cout << "current Lagrange mult for nonLin Cstr: " << lagMult_.nonLinear.transpose().format(CleanFmt) << std::endl;
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
    
    probEval_.diffLag.resize(1, problem.M().dim());

    probEval_.linearizedInfBndLinCstr.resize(cstrMngr_.totalDimLin());
    probEval_.linearizedInfBndLinCstr.resize(cstrMngr_.totalDimLin());
    probEval_.linearizedSupBndNonLinCstr.resize(cstrMngr_.totalDimNonLin());
    probEval_.linearizedSupBndNonLinCstr.resize(cstrMngr_.totalDimNonLin());

    z_.resize(problem.M().dim());
    lagMult_.linear.resize(cstrMngr_.totalDimLin());
    lagMult_.nonLinear.resize(cstrMngr_.totalDimNonLin());

    probEval_.allCstr.resize(cstrMngr_.totalDimNonLin()+cstrMngr_.totalDimLin());
    lagMult_.all.resize(cstrMngr_.totalDimNonLin()+cstrMngr_.totalDimLin());
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
      // TODO this does not seem very efficient
      p.evalLinCstr(cstrMngr_.getViewLin(probEval_.linCstr,i),i);
      p.evalLinCstrDiff(cstrMngr_.getViewLin(probEval_.diffLinCstr,i),i);
      p.getLinCstrLB(cstrMngr_.getViewLin(probEval_.linCstrLB,i),i);
      p.getLinCstrUB(cstrMngr_.getViewLin(probEval_.linCstrUB,i),i);

      p.evalNonLinCstr(cstrMngr_.getViewNonLin(probEval_.nonLinCstr,i),i);
      p.evalNonLinCstrDiff(cstrMngr_.getViewNonLin(probEval_.diffNonLinCstr,i),i);
      p.getNonLinCstrLB(cstrMngr_.getViewNonLin(probEval_.nonLinCstrLB,i),i);
      p.getNonLinCstrUB(cstrMngr_.getViewNonLin(probEval_.nonLinCstrUB,i),i);
    }
    probEval_.allCstr.head(cstrMngr_.totalDimLin()) = probEval_.linCstr;
    probEval_.allCstr.tail(cstrMngr_.totalDimNonLin()) = probEval_.nonLinCstr;

    probEval_.linearizedInfBndLinCstr = probEval_.linCstr - probEval_.linCstrLB;
    probEval_.linearizedSupBndLinCstr = probEval_.linCstr - probEval_.linCstrUB;
    probEval_.linearizedInfBndNonLinCstr = probEval_.nonLinCstr - probEval_.nonLinCstrLB;
    probEval_.linearizedSupBndNonLinCstr = probEval_.nonLinCstr - probEval_.nonLinCstrUB;

    probEval_.lag = computeLagrangian();
    probEval_.diffLag = computeDiffLagrangian();
  }

  double Solver::computeLagrangian()
  {
    double res = probEval_.obj; 
    for( Index i = 0; i<cstrMngr_.totalDimLin(); ++i)
    {
      //only the constraints that are violated appear in the lagrangian. The
      //valid ones have a Lagrange multiplier of value 0. Cf Note on
      //implementation details
      res = res + lagMult_.linear[i]*(fmin(probEval_.linearizedInfBndLinCstr[i],0)
             + fmax(probEval_.linearizedSupBndLinCstr[i], 0)); 
    }
    for( Index i = 0; i<cstrMngr_.totalDimNonLin(); ++i)
    {
      //only the constraints that are violated appear in the lagrangian. The
      //valid ones have a Lagrange multiplier of value 0. Cf Note on
      //implementation details
      res = res + lagMult_.nonLinear[i]*(fmin(probEval_.linearizedInfBndNonLinCstr[i],0)
             + fmax(probEval_.linearizedSupBndNonLinCstr[i], 0)); 
    }
    return res;
  }

  Eigen::MatrixXd Solver::computeDiffLagrangian()
  {
    Eigen::MatrixXd res(1, problem_->M().dim());
    res = probEval_.diffObj 
            + lagMult_.linear.transpose()*probEval_.diffLinCstr
            + lagMult_.nonLinear.transpose()*probEval_.diffNonLinCstr;
    return res;
  }

  bool Solver::convergence(
      double tau_P, double tau_D, const Point& x, 
      const Eigen::VectorXd& lagMult, const Eigen::VectorXd& cstr, 
      const Eigen::MatrixXd& diffLag) const
  {
    std::cout << tau_P << std::endl;
    std::cout << tau_D << std::endl;
    std::cout << x << std::endl;
    std::cout << lagMult << std::endl;
    std::cout << cstr << std::endl;
    std::cout << diffLag << std::endl;
    return true;
  }
}
