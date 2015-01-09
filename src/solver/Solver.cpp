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
    lagMult_.initOnes();
    updateAllProblemData(problem);
    printStatus();

    //int iter = 0;
    //int maxIter = 1000;
    //double epsilon_P = 1e-6;
    //double epsilon_D = 1e-6;
    bool converged = convergence(opt_.epsilon_P , opt_.epsilon_P, 
        x0, lagMult_.all, probEval_.allInfCstr, probEval_.allSupCstr, probEval_.diffLag);
    std::cout << converged << std::endl;
     

    return Results({ x0, CONVERGE, {} });
  }

  void Solver::printStatus()
  {
    std::cout << "================================================================="<< std::endl;
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    std::cout << "current x = " << problem_->x() << std::endl;
    std::cout << "current z = " << z_.transpose().format(CleanFmt) << std::endl;
    lagMult_.print();
    probEval_.print();
    std::cout << "================================================================="<< std::endl;
  }

  void Solver::initSolver(Problem& problem)
  {
    problem_ = &problem;
    opt_.maxIter = 100;
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

    probEval_.infBndCstr.resize(problem.M().dim());
    probEval_.supBndCstr.resize(problem.M().dim());
    probEval_.infLinCstr.resize(cstrMngr_.totalDimLin());
    probEval_.supLinCstr.resize(cstrMngr_.totalDimLin());
    probEval_.infNonLinCstr.resize(cstrMngr_.totalDimNonLin());
    probEval_.supNonLinCstr.resize(cstrMngr_.totalDimNonLin());
    
    probEval_.allInfCstr.resize(problem.M().dim() 
                                + cstrMngr_.totalDimNonLin() 
                                + cstrMngr_.totalDimLin());
    probEval_.allSupCstr.resize(problem.M().dim() 
                                + cstrMngr_.totalDimNonLin() 
                                + cstrMngr_.totalDimLin());
    

    z_.resize(problem.M().dim());
    lagMult_.bounds.resize(problem.M().dim());
    lagMult_.linear.resize(cstrMngr_.totalDimLin());
    lagMult_.nonLinear.resize(cstrMngr_.totalDimNonLin());

    probEval_.allCstr.resize(cstrMngr_.totalDimNonLin()+cstrMngr_.totalDimLin());
    lagMult_.all.resize(problem.M().dim() + cstrMngr_.totalDimNonLin() + cstrMngr_.totalDimLin());
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

    probEval_.infBndCstr =  z_ - probEval_.tangentLB;
    probEval_.supBndCstr =  z_ - probEval_.tangentUB;
    probEval_.infLinCstr = probEval_.linCstr - probEval_.linCstrLB;
    probEval_.supLinCstr = probEval_.linCstr - probEval_.linCstrUB;
    probEval_.infNonLinCstr = probEval_.nonLinCstr - probEval_.nonLinCstrLB;
    probEval_.supNonLinCstr = probEval_.nonLinCstr - probEval_.nonLinCstrUB;

    probEval_.allInfCstr.head(problem_->M().dim()) = probEval_.infBndCstr;
    probEval_.allInfCstr.segment(problem_->M().dim(), cstrMngr_.totalDimLin())= probEval_.infLinCstr ;
    probEval_.allInfCstr.tail(cstrMngr_.totalDimNonLin()) = probEval_.infNonLinCstr ;

    probEval_.allSupCstr.head(problem_->M().dim()) = probEval_.supBndCstr;
    probEval_.allSupCstr.segment(problem_->M().dim(), cstrMngr_.totalDimLin())= probEval_.supLinCstr ;
    probEval_.allSupCstr.tail(cstrMngr_.totalDimNonLin()) = probEval_.supNonLinCstr ;

    probEval_.allCstr.head(cstrMngr_.totalDimLin()) = probEval_.linCstr;
    probEval_.allCstr.tail(cstrMngr_.totalDimNonLin()) = probEval_.nonLinCstr;

    probEval_.lag = computeLagrangian();
    probEval_.diffLag = computeDiffLagrangian();
  }

  double Solver::computeLagrangian()
  {
    double res = probEval_.obj; 
    for( Index i = 0; i<problem_->M().dim(); ++i)
    {
      //only the constraints that are violated appear in the lagrangian. The
      //valid ones have a Lagrange multiplier of value 0. Cf Note on
      //implementation details
      res = res + lagMult_.bounds[i]*(fmin(probEval_.infBndCstr[i],0)
             + fmax(probEval_.supBndCstr[i], 0)); 
    }
    for( Index i = 0; i<cstrMngr_.totalDimLin(); ++i)
    {
      //only the constraints that are violated appear in the lagrangian. The
      //valid ones have a Lagrange multiplier of value 0. Cf Note on
      //implementation details
      res = res + lagMult_.linear[i]*(fmin(probEval_.infLinCstr[i],0)
             + fmax(probEval_.supLinCstr[i], 0)); 
    }
    for( Index i = 0; i<cstrMngr_.totalDimNonLin(); ++i)
    {
      //only the constraints that are violated appear in the lagrangian. The
      //valid ones have a Lagrange multiplier of value 0. Cf Note on
      //implementation details
      res = res + lagMult_.nonLinear[i]*(fmin(probEval_.infNonLinCstr[i],0)
             + fmax(probEval_.supNonLinCstr[i], 0)); 
    }
    return res;
  }

  Eigen::MatrixXd Solver::computeDiffLagrangian()
  {
    Eigen::MatrixXd res(1, problem_->M().dim());
    res = probEval_.diffObj 
            + lagMult_.bounds.transpose()
            + lagMult_.linear.transpose()*probEval_.diffLinCstr
            + lagMult_.nonLinear.transpose()*probEval_.diffNonLinCstr;
    return res;
  }

  bool Solver::convergence(
      double tau_P, double tau_D, const Point& x, 
      const Eigen::VectorXd& lagMult, 
      const Eigen::VectorXd& infCstr, 
      const Eigen::VectorXd& supCstr, 
      const Eigen::MatrixXd& diffLag) const
  {
    std::cout << "-----------------Convergence test---------------" << std::endl;
    std::cout << "tau_P: " << tau_P << std::endl;
    std::cout << "tau_D: " << tau_D << std::endl;
    std::cout << "x: " << x << std::endl;
    std::cout << "lagMult: " << lagMult.transpose() << std::endl;
    std::cout << "infCstr: " << infCstr.transpose() << std::endl;
    std::cout << "supCstr: " << supCstr.transpose() << std::endl;
    std::cout << "diffLag: " << diffLag << std::endl;

    Eigen::VectorXd invMapX(problem_->M().dim());
    problem_->M().invMap(invMapX, x.value());
    std::cout << "invMap(x) = " << invMapX.transpose() << std::endl;
    double tau_x = tau_P*(1+invMapX.lpNorm<Eigen::Infinity>());
    double tau_l = tau_D*(1+lagMult.lpNorm<Eigen::Infinity>());

    std::cout << "tau_x: " << tau_x << std::endl;
    std::cout << "tau_l: " << tau_l << std::endl;

    bool converged = true;
    if(!(diffLag.array().abs() <= tau_l).all())
      converged = false;
    for(Index i = 0; i<lagMult.size(); ++i)
    {
      if(!((lagMult[i]<-tau_l && fabs(infCstr(i))<tau_x)
          || (fabs(lagMult[i])>tau_l && infCstr(i)>=-tau_x && supCstr(i)<=tau_x)
          || (lagMult[i]>tau_l && fabs(supCstr(i))<tau_x)))
      {
        std::cout << "Cstr " << i << " failure" << std::endl;
        converged = false;
      }
      else
      {
        std::cout << "Cstr " << i << " success" << std::endl;
      }
    }
    std::cout << "------------------------------------------------" << std::endl;
    return converged;
  }
}
