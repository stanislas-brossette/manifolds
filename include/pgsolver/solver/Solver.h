#ifndef _PGS_SOLVER_H_
#define _PGS_SOLVER_H_

#include <iostream>
#include <Eigen/Core>

#include <pgsolver/manifolds/Point.h>

#include <pgsolver/solver/ConstraintManager.h>
#include <pgsolver/solver/Filter.h>
#include <pgsolver/solver/HessianUpdater.h>
#include <pgsolver/solver/Problem.h>
#include <pgsolver/solver/Results.h>
#include <pgsolver/solver/SolverOptions.h>
#include <pgsolver/solver/TempData.h>

#include <EigenQP/LSSOL.h>
#include <EigenQP/LSSOL_LP.h>
#include <EigenQP/LSSOL_QP.h>

namespace pgs
{
  class PGS_API Solver
  {
    public:
      Solver();
      /// \brief Solves the optimization problem described in p starting from the initial guess x0
      Results solve(Problem& p, Point& x0);
      /// \brief Displays the current evaluation of the problem that is considered by the solver
      // TODO: Add a security to this method so that it can't be called before
      // solver is initialized
      void printStatus();

    public:
      void setHessianBFGS(){opt_.hessianUpdateMethod = BFGS;};
      void setHessianSR1(){opt_.hessianUpdateMethod = SR1;};
      void setHessianEXACT(){opt_.hessianUpdateMethod = EXACT;};

      /// \brief Accessor for probEval
      const ProblemEvaluation& probEval() const { return probEval_;};

      /// \brief Updates Everything in probEval
      void updateAllProblemData(Problem& p);

      /// \brief Updates the Objective
      void updateObj(Problem& p);

      /// \brief Updates all the Constraints
      void updateAllCstr(Problem& p);

      /// \brief Computes the cstrViolation of a set of constraint values. 
      /// This method requires updated constraints
      double computeCstrViolation(
                    const Eigen::VectorXd& lb, 
                    const Eigen::VectorXd& c, 
                    const Eigen::VectorXd& ub);

      /// \brief Option for the solver
      SolverOptions opt_;

      /// \brief Computes the Cost and constraint violation given an
      /// infeasibility status vector
      Eigen::Vector2d computeFHforFilter(
          const ProblemEvaluation& p,const Eigen::VectorXi infeasStatus);

    protected:
      /// \brief Initializes the solver, makes all the memory allocations
      void initSolver(Problem& p);

      /// \brief Tests the convergence of the solver based on the criterion
      /// presented in SNOPT paper page 108
      /// \param tau_P Primal problem convergence criterion
      /// \param tau_D Dual problem convergence criterion
      /// \param x point where the convergence is tested
      /// \param lagMult Concatenation of all the Lagrange Multiplier with on
      /// top the ones associated with linear constraints, followed by the ones
      /// associated with non-linear constraints
      /// \param cstr Concatenation of all the constraints with on
      /// top the linear constraints, followed by the non-linear constraints
      /// \param diffLag Derivative of the Lagrangian (Jacobian. Should be a line-vector)
      bool convergence(
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
        const Eigen::MatrixXd& diffLag) const;
      bool KKTTestCstr(
        double tau_l, double tau_x,
        const Eigen::VectorXd& lagMult,
        const Eigen::VectorXd& infCstr,
        const Eigen::VectorXd& supCstr) const;

      /// \brief Restoration method makes sure that the quadratic approximated problem
      /// has a solution. If it doesn't, it finds a new evaluation point for
      /// which it is feasible
      void restoration();

      /// \brief Computes the cost and constraints for restoration based on
      /// feasibility results
      void computeRestorationQuantities(ProblemEvaluation& probEval, Index& nCstr, const SolverOptions& opt);

      /**
      \brief Method testing the feasibility of a set of constraints
      This method solves the following problem:
      \f{align}{
        \min \sum {{\bf v}_i} + \sum {{\bf w}_i} \nonumber
        \\ \text{subject to }&
        \left\{
        \begin{array}{lr}
          lb \leq {\bf x}\leq ub
          \\ lb \leq \nabla c_k.{\bf x} + c_k + {\bf v} - {\bf w} \leq ub
          \\ 0 \leq {\bf v}_i \leq +\infty 
          \\ 0 \leq {\bf w}_i \leq +\infty 
        \\ \end{array} \nonumber
        \right.
      \f}
      **/
      bool feasibility( ProblemEvaluation& probEval, double feasibilityMin,
                        Eigen::VectorXd& feasibleVector,
                        Eigen::VectorXd& infeasibilityInf, Eigen::VectorXd& infeasibilitySup,
                        LagrangeMultipliers& feasLagMult);

      bool feasibility( Eigen::VectorXd& xl, Eigen::VectorXd& xu,
                        Eigen::VectorXd& cvec, Eigen::MatrixXd& C,
                        Eigen::VectorXd& cl, Eigen::VectorXd& cu,
                        double feasibilityMin,
                        Eigen::VectorXd& feasibleVector,
                        Eigen::VectorXd& infeasibilityInf, Eigen::VectorXd& infeasibilitySup,
                        LagrangeMultipliers& feasLagMult);

      /// \brief Computes the value of the Lagrangian of the problem
      double computeLagrangian();
      Eigen::MatrixXd computeDiffLagrangian();

    private:
      /// \brief Structure containing The Lagrange Multiplier for Linear and nonLinear Constraints
      LagrangeMultipliers lagMult_;
      /// \brief Structure containing The Lagrange Multiplier for Linear and nonLinear Constraints
      /// that is only used for the feasibility phases
      LagrangeMultipliers feasibilityLagMult_;
      /// \brief Structure containing The Lagrange Multiplier for Linear and nonLinear Constraints
      /// that is only used for the restoration phases
      LagrangeMultipliers restorationLagMult_;
      /// \brief Vector that contains the step data
      Eigen::VectorXd z_;

      /// \brief Set of vector and matrices representing the State of the problem
      ProblemEvaluation probEval_;
      /// \brief Objet that knows and manages the memory location for all the constraints of the problem
      ConstraintManager cstrMngr_;

      /// \brief Filter
      Filter filter_;
      /// \brief Filter for restoration
      Filter restFilter_;

      /// \brief pointer on the problem considered
      Problem* problem_;

      /// \brief QP solver
      Eigen::LSSOL_QP QPSolver_;
      /// \brief LP solver
      Eigen::LSSOL_LP LPSolver_;
      /// \brief Restoration QP solver
      Eigen::LSSOL_QP RestQPSolver_;

  };
}

#endif //_PGS_SOLVER_H_
