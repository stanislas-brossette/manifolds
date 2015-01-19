#include <pgsolver/solver/ConstraintManager.h>
#include <pgsolver/solver/LineSearcher.h>
#include <pgsolver/solver/HessianUpdater.h>
#include <pgsolver/solver/Filter.h>
#include <pgsolver/solver/Solver.h>

namespace pgs
{
  Solver::Solver()
  {
    std::cout << "New Solver" << std::endl;
    filter_ = Filter();
  }

  Results Solver::solve(Problem& problem, Point& x0)
  {
    problem.setX(x0);
    cstrMngr_.init(problem);
    initSolver(problem);
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

    z_.setZero();
    lagMult_.initOnes();
    updateAllProblemData(problem);
    std::cout << "================== Initial Conditions ==========================="<< std::endl;
    printStatus();

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
      std::cout <<std::endl<< "********************Iteration " << iter <<"*********************"<< std::endl;

      //Test Feasibility
      bool feasible = feasibility(probEval_, opt_.epsilonFeasibility);
      std::cout << "Problem Feasibility = " << feasible << std::endl;

      //Resolution of the quadratic tangent problem
      QPSolver_.solve(
          probEval_.Hessian,
          probEval_.diffObj.transpose(),
          probEval_.allDiffCstr,
          static_cast<int>(probEval_.allDiffCstr.rows()),
          -probEval_.allInfCstr,
          -probEval_.allSupCstr,
          probEval_.tangentLB,
          probEval_.tangentUB);
      z_ = QPSolver_.result();

      //Globalization
      switch(opt_.globalizationMethod){
        case NONE: alpha = 1; break;
        case LINESEARCH:
          alpha = LineSearcher::LineSearch(*this, problem, filter_, z_);
          break;
        case TRUSTREGION:
          std::runtime_error("EXACT update is not implemented yet");
          break;
      }
      std::cout << "LineSearch gives alpha = " << alpha << std::endl;

      lagMult_.bounds = (1-alpha)*lagMult_.bounds + alpha*(-QPSolver_.clambda().head(lagMult_.bounds.size()));
      lagMult_.linear = (1-alpha)*lagMult_.linear + alpha*(-QPSolver_.clambda().segment(lagMult_.bounds.size(), lagMult_.linear.size()));
      lagMult_.nonLinear = (1-alpha)*lagMult_.nonLinear + alpha*(-QPSolver_.clambda().tail(lagMult_.nonLinear.size()));

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

      //printStatus();
    }
    std::cout << "=============== Solution at iteration " << iter << " ========================="<< std::endl;
    printStatus();

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

    probEval_.allInfCstr.resize(cstrMngr_.totalDimNonLin()
                                + cstrMngr_.totalDimLin());
    probEval_.allInfCstr.setZero();
    probEval_.allSupCstr.resize(cstrMngr_.totalDimNonLin()
                                + cstrMngr_.totalDimLin());
    probEval_.allSupCstr.setZero();


    z_.resize(probEval_.varDim);
    lagMult_.bounds.resize(probEval_.varDim);
    lagMult_.bounds.setZero();
    lagMult_.linear.resize(cstrMngr_.totalDimLin());
    lagMult_.linear.setZero();
    lagMult_.nonLinear.resize(cstrMngr_.totalDimNonLin());
    lagMult_.nonLinear.setZero();

    probEval_.allCstr.resize(cstrMngr_.totalDimNonLin()+cstrMngr_.totalDimLin());
    probEval_.allCstr.setZero();
    probEval_.allDiffCstr.resize(cstrMngr_.totalDim(), probEval_.varDim);
    probEval_.allDiffCstr.setZero();

    probEval_.violCstr.resize(probEval_.varDim + cstrMngr_.totalDim());
    probEval_.violCstr.setZero();

    lagMult_.all.resize(probEval_.varDim + cstrMngr_.totalDim());
    lagMult_.all.setZero();

    Index nFeasCstr = 2*cstrMngr_.totalDim();
    probEval_.feasibilityCostF.resize(probEval_.varDim + nFeasCstr );
    probEval_.feasibilityCostF.head(probEval_.varDim) = Eigen::VectorXd::Zero(probEval_.varDim);
    probEval_.feasibilityCostF.tail( nFeasCstr ) = Eigen::VectorXd::Constant( nFeasCstr, 1);
    probEval_.feasibilityAllDiffCstr.resize( nFeasCstr, probEval_.varDim + nFeasCstr );
    probEval_.feasibilityAllDiffCstr.setZero();
    probEval_.feasibilityAllDiffCstr.block(0,probEval_.varDim, nFeasCstr, nFeasCstr) = Eigen::MatrixXd::Identity(nFeasCstr, nFeasCstr);
    probEval_.feasibilityLB.resize(probEval_.varDim + nFeasCstr);
    probEval_.feasibilityLB.setZero();
    probEval_.feasibilityLB.head(probEval_.varDim) = probEval_.tangentLB;
    probEval_.feasibilityUB.resize(probEval_.varDim + nFeasCstr);
    probEval_.feasibilityUB.head(probEval_.varDim) = probEval_.tangentUB;
    probEval_.feasibilityUB.tail(nFeasCstr) = 
           Eigen::VectorXd::Constant(nFeasCstr, std::numeric_limits<double>::infinity());

    probEval_.feasibilityInfCstr.resize(nFeasCstr);
    probEval_.feasibilityInfCstr.setZero();
    probEval_.feasibilitySupCstr.resize(nFeasCstr);
    probEval_.feasibilitySupCstr.setZero();

    probEval_.infeasibilityResult.resize(nFeasCstr);
    probEval_.infeasibilityResult.setZero();

    QPSolver_ = Eigen::LSSOL(int(probEval_.varDim), int(cstrMngr_.totalDim()));

    //TODO: This is WRONG. It should be a LSSOL solver of type LP. Dunno yet how
    //to make one.
    LPSolver_ = Eigen::LSSOL(int(probEval_.varDim + nFeasCstr), int(nFeasCstr));

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
    probEval_.allCstr.head(cstrMngr_.totalDimLin()) = probEval_.linCstr;
    probEval_.allCstr.tail(cstrMngr_.totalDimNonLin()) = probEval_.nonLinCstr;

    probEval_.infLinCstr = probEval_.linCstr - probEval_.linCstrLB;
    probEval_.supLinCstr = probEval_.linCstr - probEval_.linCstrUB;
    probEval_.infNonLinCstr = probEval_.nonLinCstr - probEval_.nonLinCstrLB;
    probEval_.supNonLinCstr = probEval_.nonLinCstr - probEval_.nonLinCstrUB;

    probEval_.allInfCstr.head(cstrMngr_.totalDimLin())= probEval_.infLinCstr ;
    probEval_.allInfCstr.tail(cstrMngr_.totalDimNonLin()) = probEval_.infNonLinCstr ;

    probEval_.allSupCstr.head(cstrMngr_.totalDimLin())= probEval_.supLinCstr ;
    probEval_.allSupCstr.tail(cstrMngr_.totalDimNonLin()) = probEval_.supNonLinCstr ;

    probEval_.allCstr.head(cstrMngr_.totalDimLin()) = probEval_.linCstr;
    probEval_.allCstr.tail(cstrMngr_.totalDimNonLin()) = probEval_.nonLinCstr;

    probEval_.allDiffCstr.block(0,0, cstrMngr_.totalDimLin(), p.M().dim()) = probEval_.diffLinCstr;
    probEval_.allDiffCstr.block(cstrMngr_.totalDimLin(), 0, cstrMngr_.totalDimNonLin(), p.M().dim()) = probEval_.diffNonLinCstr;

    probEval_.lag = computeLagrangian();
    probEval_.prevDiffLag = probEval_.diffLag;
    probEval_.diffLag = computeDiffLagrangian();



    Index nFeasCstr = 2*cstrMngr_.totalDim();
    for(Index i = 0; i < cstrMngr_.totalDim(); ++i)
    {
      probEval_.feasibilityAllDiffCstr.block(2*i, 0, 1, probEval_.varDim) =
        probEval_.allDiffCstr.row(i);
      probEval_.feasibilityAllDiffCstr.block(2*i+1, 0, 1, probEval_.varDim) =
        - probEval_.allDiffCstr.row(i);
      probEval_.feasibilityInfCstr(2*i) = -probEval_.allInfCstr(i);
      probEval_.feasibilityInfCstr(2*i+1) = probEval_.allSupCstr(i);
    }
    probEval_.feasibilitySupCstr = Eigen::VectorXd::Constant(nFeasCstr, std::numeric_limits<double>::infinity());
    probEval_.feasibilityLB.head(probEval_.varDim) = probEval_.tangentLB;
    probEval_.feasibilityLB.tail(nFeasCstr) = Eigen::VectorXd::Zero(nFeasCstr);
    probEval_.feasibilityUB.head(probEval_.varDim) = probEval_.tangentUB;
    probEval_.feasibilityUB.tail(nFeasCstr) = Eigen::VectorXd::Constant(nFeasCstr, std::numeric_limits<double>::infinity());
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
  }

  void Solver::updateViolations(Problem& p)
  {
    updateAllCstr(p);

    for (Index i = 0; i < p.M().dim(); ++i)
    {
      probEval_.violCstr(i) = fmax(
          fmax(probEval_.tangentLB(i) - p.z()(i),
            p.z()(i) - probEval_.tangentUB(i)),0);
    }
    for (Index i = 0; i<probEval_.linCstr.size(); ++i)
    {
      probEval_.violCstr(i) = fmax(
          fmax(probEval_.linCstrLB(i) - probEval_.linCstr(i),
            probEval_.linCstr(i) - probEval_.linCstrUB(i)),0);
    }
    for (Index i = 0; i<probEval_.nonLinCstr.size(); ++i)
    {
      probEval_.violCstr(i) = fmax(
          fmax(probEval_.nonLinCstrLB(i) - probEval_.nonLinCstr(i),
            probEval_.nonLinCstr(i) - probEval_.nonLinCstrUB(i)),0);
    }
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
    std::cout << "Test KKT diff Lagrangian: " << std::endl;
    bool convergedLag = (diffLag.array().abs() <= tau_l).all();
    std::cout << "convergedLag = " << convergedLag << std::endl;
    //Test bounds cstr
    std::cout << "Test KKT Bound Cstr " << std::endl;
    bool convergedBounds = KKTTestCstr(tau_l, tau_x, lagMultBnd, -tangentLB, -tangentUB);
    //Test Linear cstr
    std::cout << "Test KKT Linear Cstr " << std::endl;
    bool convergedLin = KKTTestCstr(tau_l, tau_x, lagMultLin, infCstrLin, supCstrLin);
    //Test NonLinear cstr
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
        std::cout << "Cstr " << i << " failure" << std::endl;
        //std::cout << "!((lagMult[i]<-tau_l && fabs(infCstr(i))<tau_x) || (fabs(lagMult[i])>tau_l && infCstr(i)>=-tau_x && supCstr(i)<=tau_x) || (lagMult[i]>tau_l && fabs(supCstr(i))<tau_x)))"<< std::endl;
        //std::cout << "!((" << lagMult[i] << " < " << -tau_l << " && " << fabs(infCstr(i)) << " < " << tau_x << ") || ( " << fabs(lagMult[i]) << " > " << tau_l << " && " << infCstr(i)<< " >= " << -tau_x << " && " << supCstr(i) << " <= " << tau_x << ") || (" << lagMult[i] << " > " << tau_l << " && " << fabs(supCstr(i)) << " < " << tau_x << "))" << std::endl;
        converged = false;
      }
      else
      {
        std::cout << "Cstr " << i << " success" << std::endl;
      }
    }
    return converged;
  }

  bool Solver::feasibility(const ProblemEvaluation& probEval, double eps_feasibility)
  {
    //TODO Cleanup that method
    std::cout << "Feasibility Test:"<< std::endl;
    std::cout << "probEval_.feasibilityAllDiffCstr = \n" << probEval_.feasibilityAllDiffCstr << std::endl;
    std::cout << "probEval_.feasibilityInfCstr = \n" << probEval_.feasibilityInfCstr << std::endl;
    std::cout << "probEval_.feasibilitySupCstr = \n" << probEval_.feasibilitySupCstr << std::endl;
    std::cout << "probEval_.feasibilityLB = \n" << probEval_.feasibilityLB << std::endl;
    std::cout << "probEval_.feasibilityUB = \n" << probEval_.feasibilityUB << std::endl;
    Eigen::MatrixXd BigZero(probEval.feasibilityCostF.size(), probEval.feasibilityCostF.size());
    BigZero = Eigen::MatrixXd::Zero(probEval.feasibilityCostF.size(), probEval.feasibilityCostF.size());
    LPSolver_.solve(
        BigZero,
        probEval.feasibilityCostF,
        probEval.feasibilityAllDiffCstr,
        static_cast<int>(probEval.feasibilityAllDiffCstr.rows()),
        probEval.feasibilityInfCstr,
        probEval.feasibilitySupCstr,
        probEval.feasibilityLB,
        probEval.feasibilityUB);
    Index nFeasCstr = probEval.infeasibilityResult.size();
    Eigen::VectorXd LPresult(nFeasCstr);
    LPresult = LPSolver_.result();
    std::cout << "Result feasibility test:\n" << LPresult << std::endl;
    std::cout << "Result feasibility test:\n" << LPSolver_.result().tail(probEval.infeasibilityResult.size()) << std::endl;
    Eigen::VectorXd infeasibilityResult(2*probEval.linCstrDim + 2*probEval.nonLinCstrDim);
    for(Index i = 0; i < nFeasCstr; ++i)
    {
      infeasibilityResult(i) = LPresult(probEval.varDim + i);
    }
    //probEval.infeasibilityResult = result.tail(probEval.infeasibilityResult.size());    
    std::cout << "infeasibilityResult:\n" << infeasibilityResult << std::endl;
    bool result = (infeasibilityResult.lpNorm<Eigen::Infinity>() < eps_feasibility);
    return result;
  }


  const ProblemEvaluation& Solver::probEval() const
  {
    return probEval_;
  }
}
