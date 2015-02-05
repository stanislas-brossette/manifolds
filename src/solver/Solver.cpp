#include <stdexcept>
#include <limits>
#include <pgsolver/solver/ConstraintManager.h>
#include <pgsolver/solver/LineSearcher.h>
#include <pgsolver/solver/HessianUpdater.h>
#include <pgsolver/solver/Filter.h>
#include <pgsolver/solver/Solver.h>

namespace pgs
{
  Solver::Solver()
    : filter_(1e-12),
      restFilter_(1e-12)
  {
  }

  Results Solver::solve(Problem& problem, Point& x0)
  {
    problem.setX(x0);
    cstrMngr_.init(problem);
    initSolver(problem);
    z_.setZero();
    lagMult_.initOnes();
    restorationLagMult_.initOnes();
    updateAllProblemData(problem);
    //////////On init prev values are set as equal to values////////////
    probEval_.prevDiffLag = probEval_.diffLag;
    probEval_.prevDiffObj = probEval_.diffObj;
    probEval_.prevDiffNonLinCstr = probEval_.diffNonLinCstr;

    if(opt_.VERBOSE >= 1) 
    {
      std::cout << "Problem with " << cstrMngr_.totalDimLin() << " Linear cstr" << std::endl;
      std::cout << "And " << cstrMngr_.totalDimNonLin() << " NonLinear cstr" << std::endl;

      std::cout << "Hessian update method is: ";
      switch(opt_.hessianUpdateMethod){
        case BFGS: std::cout << "BFGS" << std::endl; break;
        case SR1: std::cout << "SR1" << std::endl; break;
        case EXACT: std::cout << "EXACT" << std::endl; break;
      }
      std::cout << "Hessian update Type is: ";
      switch(opt_.hessianUpdateType){
        case GROUPED: std::cout << "GROUPED" << std::endl; break;
        case INDIVIDUAL: std::cout << "INDIVIDUAL" << std::endl; break;
      }
      std::cout << "Globalization method is: ";
      switch(opt_.globalizationMethod){
        case NONE: std::cout << "NONE" << std::endl; break;
        case LINESEARCH: std::cout << "LINESEARCH" << std::endl; break;
        case TRUSTREGION: std::cout << "TRUSTREGION" << std::endl; break;
      }
      std::cout << "================== Initial Conditions ==========================="<< std::endl;
      printStatus();
    }


    int iter = 0;
    double alpha = 1;

    while(!convergence(opt_.epsilon_P , opt_.epsilon_P, problem.x(),
                        lagMult_.bounds,
                        probEval_.tangentLB,
                        probEval_.tangentUB,
                        lagMult_.linear,
                        probEval_.infLinCstr,
                        probEval_.supLinCstr,
                        lagMult_.nonLinear,
                        probEval_.infNonLinCstr,
                        probEval_.supNonLinCstr,
                        probEval_.diffLag)
                        && iter < opt_.maxIter)
    {
      iter++;
      if(opt_.VERBOSE >= 1) 
      {
        std::cout <<std::endl<< "********************Iteration " << iter <<"*********************"<< std::endl;
      }

      //std::cout << "::::::::::::::::Arguments to be sent to LSSOL :::::::::::::::" << std::endl;
      //std::cout << "probEval_.Hessian = \n" << probEval_.Hessian << std::endl;
      //std::cout << "probEval_.diffObj.transpose() = \n" << probEval_.diffObj.transpose() << std::endl;
      //std::cout << "probEval_.allDiffCstr = \n" << probEval_.allDiffCstr << std::endl;
      //std::cout << "static_cast<int>(probEval_.allDiffCstr.rows()) = \n" << static_cast<int>(probEval_.allDiffCstr.rows()) << std::endl;
      //std::cout << "probEval_.allInfCstr = \n" << probEval_.allInfCstr << std::endl;
      //std::cout << "probEval_.allSupCstr = \n" << probEval_.allSupCstr << std::endl;
      //std::cout << "probEval_.tangentLB = \n" << probEval_.tangentLB << std::endl;
      //std::cout << "probEval_.tangentUB = \n" << probEval_.tangentUB << std::endl;
      //Restoration Phase
      restoration();
              
      //Resolution of the quadratic tangent problem
      QPSolver_.solve(
          probEval_.Hessian,
          probEval_.diffObj.transpose(),
          probEval_.allDiffCstr,
          probEval_.allInfCstr,
          probEval_.allSupCstr,
          probEval_.tangentLB,
          probEval_.tangentUB);

      if (QPSolver_.inform() > 0)
      {
        QPSolver_.print_inform();
        throw std::runtime_error("QP solver FAILED!!! Damnit");
      }

      z_ = QPSolver_.result();

      //Globalization
      switch(opt_.globalizationMethod){
        case NONE: alpha = 1; break;
        case LINESEARCH:
          alpha = LineSearcher::LineSearch(*this, problem, probEval_, filter_, z_, opt_);
          break;
        case TRUSTREGION:
          throw std::runtime_error("EXACT update is not implemented yet");
          break;
      }

      lagMult_.bounds = (1-alpha)*lagMult_.bounds + alpha*(-QPSolver_.lambda().head(lagMult_.bounds.size()));
      lagMult_.linear = (1-alpha)*lagMult_.linear + alpha*(-QPSolver_.lambda().segment(lagMult_.bounds.size(), lagMult_.linear.size()));
      lagMult_.nonLinear = (1-alpha)*lagMult_.nonLinear + alpha*(-QPSolver_.lambda().tail(lagMult_.nonLinear.size()));

      //Update the value of x before update problem evaluation
      problem.setX(problem.x() + alpha*z_);
      problem.setZ(Eigen::VectorXd::Zero(z_.size()));

      //Update of the values in problemEval
      updateAllProblemData(problem);

      //Update of the Hessians
      switch(opt_.hessianUpdateType){
        case GROUPED:
          HessianUpdater::hessianUpdate(
            probEval_.Hessian, problem.x(), alpha, z_,
            probEval_.prevDiffLag, probEval_.diffLag, opt_);
          break;
        case INDIVIDUAL:
          HessianUpdater::hessianUpdateIndividually(
            probEval_.Hessian, probEval_.HessianCost, probEval_.HessiansCstr,
            lagMult_.nonLinear,
            problem.x(), alpha, z_,
            probEval_.prevDiffObj, probEval_.diffObj,
            probEval_.prevDiffNonLinCstr, probEval_.diffNonLinCstr,
            opt_);
          break;
      }

    }
    
    if(opt_.VERBOSE >= 1) 
    {
      std::cout << "=============== Solution at iteration " << iter << " ========================="<< std::endl;
      printStatus();
    }

    return Results({ x0, CONVERGE, {} });
  }

  void Solver::printStatus()
  {
    std::cout << "================================================================="<< std::endl;
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    std::cout << "current x = " << problem_->x() << std::endl;
    std::cout << "previous z = " << z_.transpose().format(CleanFmt) << "(value that was used to get to this x)"<< std::endl;
    lagMult_.print();
    probEval_.print();
    std::cout << "================================================================="<< std::endl;
  }

  void Solver::initSolver(Problem& problem)
  {
    problem_ = &problem;

    filter_.setGamma(opt_.gammaFilter);
    filter_.setOption(opt_.filterOpt);

    probEval_.varDim = problem.M().dim();
    probEval_.linCstrDim = cstrMngr_.totalDimLin();
    probEval_.nonLinCstrDim = cstrMngr_.totalDimNonLin();

    probEval_.diffObj.resize(1, probEval_.varDim);
    probEval_.prevDiffObj.resize(1, probEval_.varDim);
    probEval_.diffObj.setZero();
    probEval_.prevDiffObj.setZero();
    probEval_.tangentLB.resize(probEval_.varDim);
    probEval_.tangentLB.setZero();
    probEval_.tangentUB.resize(probEval_.varDim);
    probEval_.tangentUB.setZero();

    probEval_.linCstr.resize(cstrMngr_.totalDimLin());
    probEval_.linCstr.setZero();
    probEval_.diffLinCstr.resize(cstrMngr_.totalDimLin(), probEval_.varDim);
    probEval_.diffLinCstr.setZero();
    probEval_.linCstrLB.resize(cstrMngr_.totalDimLin());
    probEval_.linCstrLB.setZero();
    probEval_.linCstrUB.resize(cstrMngr_.totalDimLin());
    probEval_.linCstrUB.setZero();

    probEval_.nonLinCstr.resize(cstrMngr_.totalDimNonLin());
    probEval_.nonLinCstr.setZero();
    probEval_.diffNonLinCstr.resize(cstrMngr_.totalDimNonLin(), probEval_.varDim);
    probEval_.diffNonLinCstr.setZero();
    probEval_.prevDiffNonLinCstr.resize(cstrMngr_.totalDimNonLin(), probEval_.varDim);
    probEval_.prevDiffNonLinCstr.setZero();
    probEval_.nonLinCstrLB.resize(cstrMngr_.totalDimNonLin());
    probEval_.nonLinCstrLB.setZero();
    probEval_.nonLinCstrUB.resize(cstrMngr_.totalDimNonLin());
    probEval_.nonLinCstrUB.setZero();

    probEval_.Hessian.resize(probEval_.varDim, probEval_.varDim);
    probEval_.Hessian.setIdentity();

    probEval_.HessianCost.resize(probEval_.varDim, probEval_.varDim);
    probEval_.HessianCost.setIdentity();

    probEval_.HessiansCstr.resize(static_cast<size_t>(cstrMngr_.totalDimNonLin()));
    for(size_t i =0; i<static_cast<size_t>(cstrMngr_.totalDimNonLin()); ++i)
    {
      probEval_.HessiansCstr[i].resize(probEval_.varDim, probEval_.varDim);
      probEval_.HessiansCstr[i].setIdentity();
    }

    probEval_.Hessian.setIdentity();

    probEval_.diffLag.resize(1, probEval_.varDim);
    probEval_.prevDiffLag.resize(1, probEval_.varDim);
    probEval_.diffLag.setZero();
    probEval_.prevDiffLag.setZero();

    probEval_.infLinCstr.resize(cstrMngr_.totalDimLin());
    probEval_.infLinCstr.setZero();
    probEval_.supLinCstr.resize(cstrMngr_.totalDimLin());
    probEval_.supLinCstr.setZero();
    probEval_.infNonLinCstr.resize(cstrMngr_.totalDimNonLin());
    probEval_.infNonLinCstr.setZero();
    probEval_.supNonLinCstr.resize(cstrMngr_.totalDimNonLin());
    probEval_.supNonLinCstr.setZero();

    probEval_.allInfCstr.resize(cstrMngr_.totalDimNonLin() + cstrMngr_.totalDimLin());
    probEval_.allInfCstr.setZero();
    probEval_.allSupCstr.resize(cstrMngr_.totalDimNonLin() + cstrMngr_.totalDimLin());
    probEval_.allSupCstr.setZero();


    z_.resize(probEval_.varDim);
    lagMult_.bounds.resize(probEval_.varDim);
    lagMult_.bounds.setZero();
    lagMult_.linear.resize(cstrMngr_.totalDimLin());
    lagMult_.linear.setZero();
    lagMult_.nonLinear.resize(cstrMngr_.totalDimNonLin());
    lagMult_.nonLinear.setZero();
    lagMult_.all.resize(probEval_.varDim + cstrMngr_.totalDim());
    lagMult_.all.setZero();

    probEval_.allCstrUB.resize(cstrMngr_.totalDimNonLin()+cstrMngr_.totalDimLin());
    probEval_.allCstrUB.setZero();
    probEval_.allCstrLB.resize(cstrMngr_.totalDimNonLin()+cstrMngr_.totalDimLin());
    probEval_.allCstrLB.setZero();
    probEval_.allCstr.resize(cstrMngr_.totalDimNonLin()+cstrMngr_.totalDimLin());
    probEval_.allCstr.setZero();
    probEval_.allDiffCstr.resize(cstrMngr_.totalDim(), probEval_.varDim);
    probEval_.allDiffCstr.setZero();


    // -------------------- FEASIBILITY DATA ---------------------------
    probEval_.nFeasCstr = 2*cstrMngr_.totalDim();
    probEval_.feasibilityCostF.resize(probEval_.varDim + 2*cstrMngr_.totalDim() );
    probEval_.feasibilityCostF.head(probEval_.varDim) = Eigen::VectorXd::Zero(probEval_.varDim);
    probEval_.feasibilityCostF.segment(probEval_.varDim, cstrMngr_.totalDim()) = 
                            Eigen::VectorXd::Constant( cstrMngr_.totalDim(), 1);
    probEval_.feasibilityCostF.tail( cstrMngr_.totalDim() ) = 
                            Eigen::VectorXd::Constant( cstrMngr_.totalDim(), 1);
    probEval_.feasibilityAllDiffCstr.resize( cstrMngr_.totalDim(), probEval_.varDim + 2*cstrMngr_.totalDim());
    probEval_.feasibilityAllDiffCstr.setZero();
    probEval_.feasibilityAllDiffCstr.block(
        0,probEval_.varDim, cstrMngr_.totalDim(), cstrMngr_.totalDim()) = 
                  Eigen::MatrixXd::Identity(cstrMngr_.totalDim(), cstrMngr_.totalDim());
    probEval_.feasibilityAllDiffCstr.block(
        0,probEval_.varDim + cstrMngr_.totalDim(), cstrMngr_.totalDim(), cstrMngr_.totalDim()) = 
                  -Eigen::MatrixXd::Identity(cstrMngr_.totalDim(), cstrMngr_.totalDim());
    probEval_.feasibilityLB.resize(probEval_.varDim + 2*cstrMngr_.totalDim());
    probEval_.feasibilityLB.setZero();
    probEval_.feasibilityLB.head(probEval_.varDim) = probEval_.tangentLB;
    probEval_.feasibilityUB.resize(probEval_.varDim + 2*cstrMngr_.totalDim());
    probEval_.feasibilityUB.head(probEval_.varDim) = probEval_.tangentUB;
    probEval_.feasibilityUB.tail(2*cstrMngr_.totalDim()) = 
           Eigen::VectorXd::Constant(2*cstrMngr_.totalDim(), std::numeric_limits<double>::infinity());

    probEval_.feasibleValue.resize(probEval_.varDim);
    probEval_.feasibleValue.setZero();
    probEval_.infeasibilityInf.resize(cstrMngr_.totalDim());
    probEval_.infeasibilityInf.setZero();
    probEval_.infeasibilitySup.resize(cstrMngr_.totalDim());
    probEval_.infeasibilitySup.setZero();

    probEval_.infeasStatus.resize(cstrMngr_.totalDim());
    probEval_.infeasStatus.setZero();

    feasibilityLagMult_.bounds.resize(probEval_.varDim);
    feasibilityLagMult_.bounds.setZero();
    feasibilityLagMult_.linear.resize(cstrMngr_.totalDimLin());
    feasibilityLagMult_.linear.setZero();
    feasibilityLagMult_.nonLinear.resize(cstrMngr_.totalDimNonLin());
    feasibilityLagMult_.nonLinear.setZero();
    feasibilityLagMult_.all.resize(probEval_.varDim + cstrMngr_.totalDim());
    feasibilityLagMult_.all.setZero();
    // -----------------------------------------------------------------

    // -------------------- RESTORATION DATA ---------------------------
    probEval_.restorationDiffObj.resize(probEval_.varDim );
    probEval_.restorationDiffObj.setZero();
    probEval_.restorationAllDiffCstr.resize( cstrMngr_.totalDim(), probEval_.varDim);
    probEval_.restorationAllDiffCstr.setZero();
    probEval_.restorationAllInfCstr.resize(cstrMngr_.totalDim());
    probEval_.restorationAllInfCstr.setZero();
    probEval_.restorationAllSupCstr.resize(cstrMngr_.totalDim());
    probEval_.restorationAllSupCstr.setZero();
    restorationLagMult_.bounds.resize(probEval_.varDim);
    restorationLagMult_.bounds.setZero();
    restorationLagMult_.linear.resize(cstrMngr_.totalDimLin());
    restorationLagMult_.linear.setZero();
    restorationLagMult_.nonLinear.resize(cstrMngr_.totalDimNonLin());
    restorationLagMult_.nonLinear.setZero();
    restorationLagMult_.all.resize(probEval_.varDim + cstrMngr_.totalDim());
    restorationLagMult_.all.setZero();
    // -----------------------------------------------------------------

    // -------------------- SOLVERS ---------------------------
    QPSolver_.resize(int(probEval_.varDim), int(cstrMngr_.totalDim()), Eigen::lssol::eType::QP2);
    LPSolver_.resize(int(probEval_.varDim+2*cstrMngr_.totalDim()), int(cstrMngr_.totalDim()));
    // --------------------------------------------------------

  }

  void Solver::updateAllProblemData(Problem& p)
  {
    probEval_.prevDiffObj = probEval_.diffObj;
    probEval_.prevDiffNonLinCstr = probEval_.diffNonLinCstr;
    p.evalObj(probEval_.obj);
    p.evalObjDiff(probEval_.diffObj);
    p.getTangentLB(probEval_.tangentLB);
    p.getTangentUB(probEval_.tangentUB);

    for (size_t i = 0; i<p.numberOfCstr(); ++i)
    {
      // for each constraint we fill the correct lines of the matrices using
      // getView
      // TODO this does not seem very efficient. Some constant data are
      // recomputed...
      p.evalLinCstr(cstrMngr_.getViewLin(probEval_.linCstr,i),i);
      p.evalLinCstrDiff(cstrMngr_.getViewLin(probEval_.diffLinCstr,i),i);
      p.getLinCstrLB(cstrMngr_.getViewLin(probEval_.linCstrLB,i),i);
      p.getLinCstrUB(cstrMngr_.getViewLin(probEval_.linCstrUB,i),i);

      p.evalNonLinCstr(cstrMngr_.getViewNonLin(probEval_.nonLinCstr,i),i);
      p.evalNonLinCstrDiff(cstrMngr_.getViewNonLin(probEval_.diffNonLinCstr,i),i);
      p.getNonLinCstrLB(cstrMngr_.getViewNonLin(probEval_.nonLinCstrLB,i),i);
      p.getNonLinCstrUB(cstrMngr_.getViewNonLin(probEval_.nonLinCstrUB,i),i);
    }

    probEval_.infLinCstr = probEval_.linCstr - probEval_.linCstrLB;
    probEval_.supLinCstr = probEval_.linCstr - probEval_.linCstrUB;
    probEval_.infNonLinCstr = probEval_.nonLinCstr - probEval_.nonLinCstrLB;
    probEval_.supNonLinCstr = probEval_.nonLinCstr - probEval_.nonLinCstrUB;

    probEval_.allInfCstr.head(cstrMngr_.totalDimLin())= -probEval_.infLinCstr ;
    probEval_.allInfCstr.tail(cstrMngr_.totalDimNonLin()) = -probEval_.infNonLinCstr ;

    probEval_.allSupCstr.head(cstrMngr_.totalDimLin())= -probEval_.supLinCstr ;
    probEval_.allSupCstr.tail(cstrMngr_.totalDimNonLin()) = -probEval_.supNonLinCstr ;

    probEval_.allCstr.head(cstrMngr_.totalDimLin()) = probEval_.linCstr;
    probEval_.allCstr.tail(cstrMngr_.totalDimNonLin()) = probEval_.nonLinCstr;
    probEval_.allCstrLB.head(probEval_.linCstrLB.size()) = probEval_.linCstrLB;
    probEval_.allCstrLB.tail(probEval_.nonLinCstrLB.size()) = probEval_.nonLinCstrLB;
    probEval_.allCstrUB.head(probEval_.linCstrUB.size()) = probEval_.linCstrUB;
    probEval_.allCstrUB.tail(probEval_.nonLinCstrUB.size()) = probEval_.nonLinCstrUB;

    probEval_.allDiffCstr.block(0,0, cstrMngr_.totalDimLin(), p.M().dim()) = probEval_.diffLinCstr;
    probEval_.allDiffCstr.block(cstrMngr_.totalDimLin(), 0, cstrMngr_.totalDimNonLin(), p.M().dim()) = probEval_.diffNonLinCstr;

    probEval_.lag = computeLagrangian();
    probEval_.prevDiffLag = probEval_.diffLag;
    probEval_.diffLag = computeDiffLagrangian();


    //TODO: Avoid those 3 lines by using a map on those matrix. They are just
    //copies of others...
    probEval_.feasibilityAllDiffCstr.block(0, 0, cstrMngr_.totalDim(), probEval_.varDim) = probEval_.allDiffCstr;
    probEval_.feasibilityLB.head(probEval_.varDim) = probEval_.tangentLB;
    probEval_.feasibilityUB.head(probEval_.varDim) = probEval_.tangentUB;

    for (Index i = 0; i < cstrMngr_.totalDim(); ++i)
    {
      assert ( !(probEval_.allInfCstr(i) == std::numeric_limits<double>::infinity() 
          && probEval_.allSupCstr(i) == -std::numeric_limits<double>::infinity()) 
          && "Both bounds of this linear constraint are -infinity and +infinity. ALWAYS feasible");
    }
  }

  void Solver::updateObj(Problem& p)
  {
    p.evalObj(probEval_.obj);
  }

  void Solver::updateAllCstr(Problem& p)
  {
    p.getTangentLB(probEval_.tangentLB);
    p.getTangentUB(probEval_.tangentUB);

    for (size_t i = 0; i<p.numberOfCstr(); ++i)
    {
      // for each constraint we fill the correct lines of the matrices using
      // getView
      // TODO this does not seem very efficient
      p.evalLinCstr(cstrMngr_.getViewLin(probEval_.linCstr,i),i);
      p.getLinCstrLB(cstrMngr_.getViewLin(probEval_.linCstrLB,i),i);
      p.getLinCstrUB(cstrMngr_.getViewLin(probEval_.linCstrUB,i),i);

      p.evalNonLinCstr(cstrMngr_.getViewNonLin(probEval_.nonLinCstr,i),i);
      p.getNonLinCstrLB(cstrMngr_.getViewNonLin(probEval_.nonLinCstrLB,i),i);
      p.getNonLinCstrUB(cstrMngr_.getViewNonLin(probEval_.nonLinCstrUB,i),i);
    }
    probEval_.allCstr.head(probEval_.linCstr.size()) = probEval_.linCstr;
    probEval_.allCstr.tail(probEval_.nonLinCstr.size()) = probEval_.nonLinCstr;
    probEval_.allCstrLB.head(probEval_.linCstrLB.size()) = probEval_.linCstrLB;
    probEval_.allCstrLB.tail(probEval_.nonLinCstrLB.size()) = probEval_.nonLinCstrLB;
    probEval_.allCstrUB.head(probEval_.linCstrUB.size()) = probEval_.linCstrUB;
    probEval_.allCstrUB.tail(probEval_.nonLinCstrUB.size()) = probEval_.nonLinCstrUB;
  }

  double Solver::computeCstrViolation(
      const Eigen::VectorXd& lb, const Eigen::VectorXd& c, const Eigen::VectorXd& ub)
  {
    assert(lb.size() == ub.size());
    assert(c.size() == ub.size());
    double violation = 0;
    for (Index i = 0; i < c.size(); ++i)
    {
      violation = fmax(fmax(lb(i) - c(i), c(i) - ub(i)),0);
    }
    return violation;
  }

  double Solver::computeLagrangian()
  {
    double res = probEval_.obj;
    for( Index i = 0; i<problem_->M().dim(); ++i)
    {
      //only the constraints that are violated appear in the lagrangian. The
      //valid ones have a Lagrange multiplier of value 0. Cf Note on
      //implementation details
      res = res + lagMult_.bounds[i]*(fmin(-probEval_.tangentLB[i],0)
             + fmax(-probEval_.tangentUB[i], 0));
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
      res = res + lagMult_.nonLinear[i]*(fmin(probEval_.infNonLinCstr[i],0.0)
             + fmax(probEval_.supNonLinCstr[i], 0.0));
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
      const Eigen::VectorXd& lagMultBnd,
      const Eigen::VectorXd& tangentLB,
      const Eigen::VectorXd& tangentUB,
      const Eigen::VectorXd& lagMultLin,
      const Eigen::VectorXd& infCstrLin,
      const Eigen::VectorXd& supCstrLin,
      const Eigen::VectorXd& lagMultNonLin,
      const Eigen::VectorXd& infCstrNonLin,
      const Eigen::VectorXd& supCstrNonLin,
      const Eigen::MatrixXd& diffLag) const
  {
    //std::cout << "-----------------Convergence test---------------" << std::endl;
    //std::cout << "tau_P: " << tau_P << std::endl;
    //std::cout << "tau_D: " << tau_D << std::endl;
    //std::cout << "x: " << x << std::endl;
    //std::cout << "lagMultBnd: " << lagMultBnd.transpose() << std::endl;
    //std::cout << "tangentLB: " << tangentLB.transpose() << std::endl;
    //std::cout << "tangentUB: " << tangentUB.transpose() << std::endl;
    //std::cout << "lagMultLin: " << lagMultLin.transpose() << std::endl;
    //std::cout << "infCstrLin: " << infCstrLin.transpose() << std::endl;
    //std::cout << "supCstrLin: " << supCstrLin.transpose() << std::endl;
    //std::cout << "lagMultNonLin: " << lagMultNonLin.transpose() << std::endl;
    //std::cout << "infCstrNonLin: " << infCstrNonLin.transpose() << std::endl;
    //std::cout << "supCstrNonLin: " << supCstrNonLin.transpose() << std::endl;
    //std::cout << "diffLag: " << diffLag << std::endl;

    Eigen::VectorXd invMapX(x.getManifold().dim());
    x.getManifold().invMap(invMapX, x.value());
    //std::cout << "invMap(x) = " << invMapX.transpose() << std::endl;
    double tau_x = tau_P*(1+invMapX.lpNorm<Eigen::Infinity>());
    double tau_l = tau_D*(1+fmax(fmax(lagMultBnd.lpNorm<Eigen::Infinity>(),lagMultLin.lpNorm<Eigen::Infinity>()),lagMultNonLin.lpNorm<Eigen::Infinity>()));

    //std::cout << "tau_x: " << tau_x << std::endl;
    //std::cout << "tau_l: " << tau_l << std::endl;

    //Test Lagrangian's derivative
    if(opt_.VERBOSE >= 1)
      std::cout << "Test KKT diff Lagrangian: " << std::endl;
    bool convergedLag = (diffLag.array().abs() <= tau_l).all();
    if(opt_.VERBOSE >= 1)
      std::cout << "convergedLag = " << convergedLag << std::endl;
    //Test bounds cstr
    if(opt_.VERBOSE >= 1)
      std::cout << "Test KKT Bound Cstr " << std::endl;
    bool convergedBounds = KKTTestCstr(tau_l, tau_x, lagMultBnd, -tangentLB, -tangentUB);
    //Test Linear cstr
    if(opt_.VERBOSE >= 1)
      std::cout << "Test KKT Linear Cstr " << std::endl;
    bool convergedLin = KKTTestCstr(tau_l, tau_x, lagMultLin, infCstrLin, supCstrLin);
    //Test NonLinear cstr
    if(opt_.VERBOSE >= 1)
      std::cout << "Test KKT NonLinear Cstr " << std::endl;
    bool convergedNonLin = KKTTestCstr(tau_l, tau_x, lagMultNonLin, infCstrNonLin, supCstrNonLin);

    bool converged = convergedLag && convergedBounds && convergedLin && convergedNonLin;

    //std::cout << "------------------------------------------------" << std::endl;
    return converged;
  }

  bool Solver::KKTTestCstr(
      double tau_l, double tau_x,
      const Eigen::VectorXd& lagMult,
      const Eigen::VectorXd& infCstr,
      const Eigen::VectorXd& supCstr) const
  {
    bool converged = true;
    for(Index i = 0; i<lagMult.size(); ++i)
    {
      if(!((lagMult[i]<-tau_l && fabs(infCstr(i))<tau_x)
          || (fabs(lagMult[i])<=tau_l && infCstr(i)>=-tau_x && supCstr(i)<=tau_x)
          || (lagMult[i]>tau_l && fabs(supCstr(i))<tau_x)))
      {
        if(opt_.VERBOSE >= 1)
          std::cout << "Cstr " << i << " failure" << std::endl;
        //std::cout << "!((lagMult[i]<-tau_l && fabs(infCstr(i))<tau_x) || (fabs(lagMult[i])>tau_l && infCstr(i)>=-tau_x && supCstr(i)<=tau_x) || (lagMult[i]>tau_l && fabs(supCstr(i))<tau_x)))"<< std::endl;
        //std::cout << "!((" << lagMult[i] << " < " << -tau_l << " && " << fabs(infCstr(i)) << " < " << tau_x << ") || ( " << fabs(lagMult[i]) << " > " << tau_l << " && " << infCstr(i)<< " >= " << -tau_x << " && " << supCstr(i) << " <= " << tau_x << ") || (" << lagMult[i] << " > " << tau_l << " && " << fabs(supCstr(i)) << " < " << tau_x << "))" << std::endl;
        converged = false;
      }
      else
      {
        if(opt_.VERBOSE >= 1)
          std::cout << "Cstr " << i << " success" << std::endl;
      }
    }
    return converged;
  }

  //bool Solver::feasibility( ProblemEvaluation& probEval, double eps_feasibility, 
  //    Eigen::VectorXd& feasibleVector, 
  //    Eigen::VectorXd& infeasibilityInf, Eigen::VectorXd& infeasibilitySup,
  //    LagrangeMultipliers& feasLagMult)
  //{
  //  if(opt_.VERBOSE >= 2)
  //    std::cout << "---------------- Feasibility -----------------"<< std::endl;

  //  LPSolver_.reset();
  //  LPSolver_.solve(
  //      probEval.feasibilityLB,
  //      probEval.feasibilityUB,
  //      probEval.feasibilityCostF,
  //      probEval.feasibilityAllDiffCstr,
  //      probEval_.allInfCstr,
  //      probEval_.allSupCstr
  //      );
  //  LPSolver_.print_result();
  //  LPSolver_.print_inform();
  //  if (!(LPSolver_.inform() == 0 || LPSolver_.inform() == 1))
  //  {
  //    LPSolver_.print_inform();
  //    std::cout << "LSSOL.istate()" << LPSolver_.istate() << std::endl;
  //    std::cout << "========Testing same problem with new solver========" << std::endl;
  //    Eigen::LSSOL_LP testLP(7,2);
  //    testLP.printLevel(1);
  //    testLP.solve(
  //      probEval.feasibilityLB,
  //      probEval.feasibilityUB,
  //      probEval.feasibilityCostF,
  //      probEval.feasibilityAllDiffCstr,
  //      probEval_.allInfCstr,
  //      probEval_.allSupCstr
  //      );
  //    testLP.print_inform();
  //    std::cout << "====================================================" << std::endl;

  //    throw std::runtime_error("Feasibility LP solver FAILED!!! Damnit");
  //  }
  //  feasibleVector = LPSolver_.result().head(probEval.varDim);
  //  infeasibilityInf = LPSolver_.result().segment(probEval.varDim, cstrMngr_.totalDim());
  //  infeasibilitySup = LPSolver_.result().tail(cstrMngr_.totalDim());
  //  bool result = (infeasibilityInf.lpNorm<Eigen::Infinity>() <= eps_feasibility) && (infeasibilitySup.lpNorm<Eigen::Infinity>() <= eps_feasibility);

  //  feasLagMult.bounds = -LPSolver_.lambda().head(feasLagMult.bounds.size());
  //  feasLagMult.linear = -LPSolver_.lambda().segment(feasLagMult.bounds.size(), feasLagMult.linear.size());
  //  feasLagMult.nonLinear = -LPSolver_.lambda().tail(feasLagMult.nonLinear.size());


  //  if(opt_.VERBOSE >= 2)
  //  {
  //    std::cout << "probEval_.feasibilityAllDiffCstr = \n" << 
  //                        probEval_.feasibilityAllDiffCstr << std::endl;
  //    std::cout << "probEval_.allInfCstr = \n" << probEval_.allInfCstr << std::endl;
  //    std::cout << "probEval_.allSupCstr = \n" << probEval_.allSupCstr << std::endl;
  //    std::cout << "probEval_.feasibilityLB = \n" << probEval_.feasibilityLB << std::endl;
  //    std::cout << "probEval_.feasibilityUB = \n" << probEval_.feasibilityUB << std::endl;
  //    std::cout << "LPSolver.result = \n" << LPSolver_.result() << std::endl;
  //    std::cout << "feasibleVector = \n" << feasibleVector << std::endl;
  //    std::cout << "infeasibilityInf = \n" << infeasibilityInf << std::endl;
  //    std::cout << "infeasibilitySup = \n" << infeasibilitySup << std::endl;
  //    std::cout << "LPSolver_.lambda() = \n" << LPSolver_.lambda() << std::endl;
  //    std::cout << "----------------------------------------------"<< std::endl;
  //  }
  //  return result;
  //}

  bool Solver::feasibility( Eigen::VectorXd& xl, Eigen::VectorXd& xu,
                            Eigen::VectorXd& cvec, Eigen::MatrixXd& C,
                            Eigen::VectorXd& cl, Eigen::VectorXd& cu,
                            double eps_feasibility,
                            Eigen::VectorXd& feasibleVector,
                            Eigen::VectorXd& infeasibilityInf, Eigen::VectorXd& infeasibilitySup,
                            LagrangeMultipliers& feasLagMult)
  {
    if(opt_.VERBOSE >= 2)
      std::cout << "---------------- Feasibility -----------------"<< std::endl;

    LPSolver_.printLevel(1004);
    LPSolver_.solve(xl, xu, cvec, C, cl, cu);
    LPSolver_.print_inform();
    LPSolver_.print_result();

    if (!(LPSolver_.inform() == 0 || LPSolver_.inform() == 1))
    {
      LPSolver_.inform();
      std::cout << "LSSOL.istate()" << LPSolver_.istate().transpose() << std::endl;
      LPSolver_.print_inform();
      std::cout << "LSSOL.istate()" << LPSolver_.istate() << std::endl;
      std::cout << "========Testing same problem with new solver========" << std::endl;
      Eigen::LSSOL_LP testLP(7,2);
      testLP.printLevel(1003);
      testLP.solve(xl, xu, cvec, C, cl, cu);
      testLP.print_inform();
      std::cout << "====================================================" << std::endl;
      throw std::runtime_error("Feasibility LP solver FAILED!!! Damnit");
    }
    feasibleVector = LPSolver_.result().head(probEval_.varDim);
    infeasibilityInf = LPSolver_.result().segment(probEval_.varDim, cstrMngr_.totalDim());
    infeasibilitySup = LPSolver_.result().tail(cstrMngr_.totalDim());
    bool result = (infeasibilityInf.lpNorm<Eigen::Infinity>() <= eps_feasibility) && (infeasibilitySup.lpNorm<Eigen::Infinity>() <= eps_feasibility);

    feasLagMult.bounds = -LPSolver_.lambda().head(feasLagMult.bounds.size());
    feasLagMult.linear = -LPSolver_.lambda().segment(feasLagMult.bounds.size(), feasLagMult.linear.size());
    feasLagMult.nonLinear = -LPSolver_.lambda().tail(feasLagMult.nonLinear.size());


    if(opt_.VERBOSE >= 2)
    {
      std::cout << "probEval_.feasibilityAllDiffCstr = \n" << 
                          probEval_.feasibilityAllDiffCstr << std::endl;
      std::cout << "probEval_.allInfCstr = \n" << probEval_.allInfCstr.transpose() << std::endl;
      std::cout << "probEval_.allSupCstr = \n" << probEval_.allSupCstr.transpose() << std::endl;
      std::cout << "probEval_.feasibilityLB = \n" << probEval_.feasibilityLB.transpose() << std::endl;
      std::cout << "probEval_.feasibilityUB = \n" << probEval_.feasibilityUB.transpose() << std::endl;
      std::cout << "LPSolver.result = \n" << LPSolver_.result().transpose() << std::endl;
      std::cout << "feasibleVector = \n" << feasibleVector.transpose() << std::endl;
      std::cout << "infeasibilityInf = \n" << infeasibilityInf.transpose() << std::endl;
      std::cout << "infeasibilitySup = \n" << infeasibilitySup.transpose() << std::endl;
      std::cout << "LPSolver_.lambda() = \n" << LPSolver_.lambda().transpose() << std::endl;
      std::cout << "----------------------------------------------"<< std::endl;
    }
    return result;
  }

  void Solver::restoration()
  {
    //TODO: In restoration, the lagMults.bounds and .linear are useless it
    //seems. 
    std::cout << "########## Restoration Phase ############ " << std::endl;

    //New filter for restoration phase
    restFilter_.reset();

    //Test Feasibility
    bool feasible = false; 
    double alphaRest = 1;

    //Reminder: feasibility puts the lagrange multipliers directly in the
    //restorationLagMult vector without ponderation by alpha for initialization
    feasible = feasibility( probEval_.feasibilityLB, probEval_.feasibilityUB,
        probEval_.feasibilityCostF, probEval_.feasibilityAllDiffCstr,
        probEval_.allInfCstr, probEval_.allSupCstr, 
        opt_.epsilonFeasibility,  probEval_.feasibleValue, 
        probEval_.infeasibilityInf, probEval_.infeasibilitySup,
        restorationLagMult_);

    if(opt_.VERBOSE >= 1) 
    {
      if (feasible)
        std::cout << "Problem Feasible. No need for restoration" << std::endl;
      else
        std::cout << "Problem NOT Feasible. NEED restoration" << std::endl;
    }
    
    int iterRest = 0;
    int maxIterRest = 100;
    Index indexCstr = 0;

    if (!feasible)
    {
      computeRestorationQuantities(probEval_, indexCstr, opt_);
      //Add the initial point to the filter
      restFilter_.add(computeFHforFilter(probEval_, probEval_.infeasStatus));

      HessianUpdater::hessianUpdateIndividually(
        probEval_.Hessian, probEval_.HessianCost, probEval_.HessiansCstr,
        restorationLagMult_.nonLinear,
        problem_->x(), alphaRest, z_,
        probEval_.prevDiffObj, probEval_.diffObj,
        probEval_.prevDiffNonLinCstr, probEval_.diffNonLinCstr,
        opt_);
    }
    RestQPSolver_.printLevel(0);
    while (!feasible && iterRest < maxIterRest)
    {
      iterRest += 1;

      RestQPSolver_.resize(int(probEval_.varDim), int(indexCstr), Eigen::lssol::QP2);

      if(opt_.VERBOSE >= 2)
      {
        std::cout << "iterRest = " << iterRest << std::endl;
        std::cout << "@@@@@@@@ Args for Restoration solver @@@@@@@@@@@"<< std::endl;
        std::cout << "probEval_.Hessian = \n" << probEval_.Hessian << std::endl;
        std::cout << "probEval_.restorationDiffObj = \n" << probEval_.restorationDiffObj << std::endl;
        std::cout << "probEval_.restorationAllDiffCstr.topRows(indexCstr) = \n" << probEval_.restorationAllDiffCstr.topRows(indexCstr) << std::endl;
        std::cout << "probEval_.restorationAllInfCstr.topRows(indexCstr) = \n" << probEval_.restorationAllInfCstr.topRows(indexCstr) << std::endl;
        std::cout << "probEval_.restorationAllSupCstr.topRows(indexCstr) = \n" << probEval_.restorationAllSupCstr.topRows(indexCstr) << std::endl;
        std::cout << "probEval_.tangentLB = \n" << probEval_.tangentLB << std::endl;
        std::cout << "probEval_.tangentUB = \n" << probEval_.tangentUB << std::endl;
        std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<< std::endl;
      }

      //TODO Hessian needs to be H = Sum(H_unFeas) + Sum(Lambda*H_feas)
      RestQPSolver_.solve(
          probEval_.tangentLB,
          probEval_.tangentUB,
          probEval_.Hessian,
          probEval_.restorationDiffObj,
          probEval_.restorationAllDiffCstr.topRows(indexCstr),
          probEval_.restorationAllInfCstr.topRows(indexCstr),
          probEval_.restorationAllSupCstr.topRows(indexCstr)
          );

      if (!(RestQPSolver_.inform() == 0 || RestQPSolver_.inform() == 1))
      {
        RestQPSolver_.inform();
        std::cout << "LSSOL.istate()" << RestQPSolver_.istate() << std::endl;
        throw std::runtime_error("Restoration QP solver FAILED!!! Damnit");
      }

      z_ = RestQPSolver_.result();
      alphaRest = LineSearcher::LineSearch(*this, *problem_, probEval_, 
                                restFilter_, z_, opt_, probEval_.infeasStatus);

      if (opt_.VERBOSE >= 1)
      {
        std::cout << "--- Restoration step " << iterRest << " ---" << std::endl;
        std::cout << "alpha*step = \n" << alphaRest*z_.transpose() << std::endl;
      }
      problem_->setX(problem_->x() + alphaRest*z_);
      problem_->setZ(Eigen::VectorXd::Zero(z_.size()));

      updateAllProblemData(*problem_);

      feasibility( probEval_.feasibilityLB, probEval_.feasibilityUB,
        probEval_.feasibilityCostF, probEval_.feasibilityAllDiffCstr,
        probEval_.allInfCstr, probEval_.allSupCstr, 
        opt_.epsilonFeasibility,  probEval_.feasibleValue, 
        probEval_.infeasibilityInf, probEval_.infeasibilitySup,
        restorationLagMult_);

      if (feasible)
      {
        std::cout << "######################################### " << std::endl;
        return;
      }

      //Approximated problem is not feasible. Restoration phase
      indexCstr = 0;
      computeRestorationQuantities(probEval_, indexCstr, opt_);

      HessianUpdater::hessianUpdateIndividually(
        probEval_.Hessian, probEval_.HessianCost, probEval_.HessiansCstr,
        restorationLagMult_.nonLinear,
        problem_->x(), alphaRest, z_,
        probEval_.prevDiffObj, probEval_.diffObj,
        probEval_.prevDiffNonLinCstr, probEval_.diffNonLinCstr,
        opt_);
      restorationLagMult_.bounds = (1-alphaRest)*restorationLagMult_.bounds + 
            alphaRest*feasibilityLagMult_.bounds;
      restorationLagMult_.linear = (1-alphaRest)*restorationLagMult_.linear + 
            alphaRest*feasibilityLagMult_.linear;
      restorationLagMult_.nonLinear = (1-alphaRest)*restorationLagMult_.nonLinear + 
            alphaRest*feasibilityLagMult_.nonLinear;
    }
    std::cout << "############ END OF RESTORATION ######### " << std::endl;
  }

  Eigen::Vector2d Solver::computeFHforFilter(
          const ProblemEvaluation& probE,const Eigen::VectorXi infeasStatus)
  {
    //TODO: I have a doubt about the value of z_ here
    double F = 0.0;
    double H = 0.0;
    Eigen::Vector2d FH;
    if (infeasStatus.size() == 0)
    {
      updateObj(*problem_);
      F = probE.obj;
      H = computeCstrViolation(probE.tangentLB, z_, probE.tangentUB);
      H += computeCstrViolation(probE.linCstrLB, probE.linCstr, probE.linCstrUB);
      H += computeCstrViolation(probE.nonLinCstrLB, probE.nonLinCstr, probE.nonLinCstrUB);
    }
    else
    {
      F = 0;
      H = computeCstrViolation(probE.tangentLB, z_, probE.tangentUB);
      for (Index i = 0; i < infeasStatus.size(); ++i)
      {
        switch (infeasStatus(i)){
          case VIOLATED_LB:
            F += fmax(0, probE.allCstrLB(i) - probE.allCstr(i));
            H += fmax(0, probE.allCstr(i) - probE.allCstrUB(i));
            break;
          case VIOLATED_UB:
            H += fmax(0, probE.allCstrLB(i) - probE.allCstr(i));
            F += fmax(0, probE.allCstr(i) - probE.allCstrUB(i));
            break;
          case SATISFIED:
            H += fmax(0, probE.allCstrLB(i) - probE.allCstr(i));
            H += fmax(0, probE.allCstr(i) - probE.allCstrUB(i));
            break;
        }
      }
    }
    FH << F, H;
    return FH;
  }

  void Solver::computeRestorationQuantities(ProblemEvaluation& probEval, Index& nCstr, 
      const SolverOptions& opt)
  {
    for (Index i = 0; i < cstrMngr_.totalDim(); ++i)
    {
      if(probEval.infeasibilityInf(i) <= opt.epsilonFeasibility &&
          probEval.infeasibilitySup(i) <= opt.epsilonFeasibility)
      {
        //Constraint i is feasible on both sides. It stays
        probEval.infeasStatus(i) = SATISFIED;
        probEval.restorationAllDiffCstr.row(nCstr) = probEval.allDiffCstr.row(i);
        probEval.restorationAllInfCstr(nCstr) = probEval.allInfCstr(i);
        probEval.restorationAllSupCstr(nCstr) = probEval.allSupCstr(i);
        nCstr += 1;
      }
      else if(probEval.infeasibilityInf(i) > opt.epsilonFeasibility &&
          probEval.infeasibilitySup(i) <= opt.epsilonFeasibility)
      {
        //Constraint i Inferior bound is not feasible (it goes to cost)
        //Constraint i Superior bound is feasible (it stays, except if infinity)
        probEval.infeasStatus(i) = VIOLATED_LB;
        probEval.restorationDiffObj -= probEval.allDiffCstr.row(i);

        if (-probEval.allSupCstr(i) != std::numeric_limits<double>::infinity())
        {
          probEval.restorationAllDiffCstr.row(nCstr) = probEval.allDiffCstr.row(i);
          probEval.restorationAllInfCstr(nCstr) = -std::numeric_limits<double>::infinity();
          probEval.restorationAllSupCstr(nCstr) = probEval.allSupCstr(i);
          nCstr += 1;
        }
      }
      else if(probEval.infeasibilityInf(i) <= opt.epsilonFeasibility &&
          probEval.infeasibilitySup(i) > opt.epsilonFeasibility)
      {
        //Constraint i Inferior bound is feasible (it stays, except if -infinity)
        //Constraint i Superior bound is not feasible (it goes to cost)
        probEval.infeasStatus(i) = VIOLATED_UB;
        probEval.restorationDiffObj += probEval.allDiffCstr.row(i);

        if (probEval.allInfCstr(i) != std::numeric_limits<double>::infinity())
        {
          probEval.restorationAllDiffCstr.row(nCstr) = probEval.allDiffCstr.row(i);
          probEval.restorationAllInfCstr(nCstr) = probEval.allInfCstr(i);
          probEval.restorationAllSupCstr(nCstr) = std::numeric_limits<double>::infinity();
          nCstr += 1;
        }
      }
      else
      {
        std::cout << "Constraint " << i << std::endl;
        std::cout << "infeasibilityInf(i) = \n" << probEval.infeasibilityInf(i) << std::endl;
        std::cout << "infeasibilitySup(i) = \n" << probEval.infeasibilitySup(i) << std::endl;
        std::cout << "allDiffCstr.row(i) = \n" << probEval.allDiffCstr.row(i) << std::endl;
        std::cout << "allInfCstr(i) = \n" << probEval.allInfCstr(i) << std::endl;
        std::cout << "allSupCstr(i) = \n" << probEval.allSupCstr(i) << std::endl;

        assert(0 && "This situation should not have happened");
      }
    }
  }
}
