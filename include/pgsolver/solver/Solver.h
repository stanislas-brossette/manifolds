#ifndef _PGS_SOLVER_H_
#define _PGS_SOLVER_H_

#include <iostream>
#include <Eigen/Core>

#include <pgsolver/solver/Problem.h>
#include <pgsolver/solver/Results.h>
#include <pgsolver/manifolds/Point.h>
#include <pgsolver/solver/TempData.h>
#include <pgsolver/solver/ConstraintManager.h>

namespace pgs
{
  class Solver
  {
    public:
      Solver();
      /// \brief Solves the optimization problem described in p starting from the initial guess x0
      Results solve(Problem& p, Point& x0);
      /// \brief Displays the current evaluation of the problem that is considered by the solver
      // TODO: Add a security to this method so that it can't be called before
      // solver is initialized
      void printStatus();

    protected:
      /// \brief Initializes the solver, makes all the memory allocations
      void initSolver(Problem& p); 
      void updateAllProblemData(Problem& p);
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
        const Eigen::VectorXd& lagMult, 
        const Eigen::VectorXd& infCstr, 
        const Eigen::VectorXd& supCstr, 
        const Eigen::MatrixXd& diffLag) const;

      /// \brief Computes the value of the Lagrangian of the problem
      double computeLagrangian();
      Eigen::MatrixXd computeDiffLagrangian();

    private:
      /// \brief Structure containing The Lagrange Multiplier for Linear and nonLinear Constraints
      LagrangeMultipliers lagMult_;    
      /// \brief Vector that contains the step data
      Eigen::VectorXd z_;             

      /// \brief Set of vector and matrices representing the State of the problem
      ProblemEvaluation probEval_;    
      /// \brief Option for the solver
      SolverOption opt_;              
      /// \brief Objet that knows and manages the memory location for all the constraints of the problem
      ConstraintManager cstrMngr_;

      /// \brief pointer on the problem considered
      Problem* problem_;
  };
}

#endif //_PGS_SOLVER_H_
