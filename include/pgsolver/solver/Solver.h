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
      void printStatus();

    private:
      /// \brief Initializes the solver, makes all the memory allocations
      void initSolver(Problem& p); 
      void updateAllProblemData(Problem& p);

    private:
      /// \brief The Lagrange Multiplier for Linear Constraints
      Eigen::VectorXd lagMultLin_;    
      /// \brief The Lagrange Multiplier for Non Linear Constraints
      Eigen::VectorXd lagMultNonLin_; 
      /// \brief Vector that contains the step data
      Eigen::VectorXd z_;             

      /// \brief Set of vector and matrices representing the State of the problem
      ProblemEvaluation probEval_;    
      /// \brief Option for the solver
      SolverOption opt_;              
      /// \brief Objet that knows and manages the memory location for all the constraints of the problem
      ConstraintManager cstrMngr_;    
  };
}

#endif //_PGS_SOLVER_H_
