#ifndef _PGS_TEMPDATA_H_
#define _PGS_TEMPDATA_H_

#include <Eigen/Core>

namespace pgs
{
  struct ProblemEvaluation
  {
    double obj;                           // Value of the costFunction
    Eigen::MatrixXd diffObj;              // Value ot the differential of the costFunction  
    Eigen::VectorXd tangentLB;            // Tangent Lower Bound
    Eigen::VectorXd tangentUB;            // Tangent Upper Bound

    Eigen::VectorXd linCstr;              // Value of the Linear Constraints
    Eigen::MatrixXd diffLinCstr;          // Value of the differential of the Linear constraints
    Eigen::VectorXd linCstrLB;            // Lower Bounds on the Linear Constraints 
    Eigen::VectorXd linCstrUB;            // Upper Bounds on the Linear Constraints

    Eigen::VectorXd nonLinCstr;           // Value of the NonLinear Constraints
    Eigen::MatrixXd diffNonLinCstr;       // Value of the differential of the NonLinear constraints
    Eigen::VectorXd nonLinCstrLB;         // Lower Bounds on the NonLinear Constraints 
    Eigen::VectorXd nonLinCstrUB;         // Upper Bounds on the NonLinear Constraints
    
    Eigen::MatrixXd Hessian;              // Hessian of the Lagrangian

    void print()
    {
      Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
      std::cout << "Current status of the problem evaluation" << std::endl;

      std::cout << std::endl << "Objective Function:" << std::endl;
      std::cout << "f(phi_x(z))=" << obj << std::endl;
      std::cout << "grad_z(f(phi_x(z))=" << diffObj.transpose().format(CleanFmt) << std::endl;

      std::cout << std::endl << "Bounds:" << std::endl;
      std::cout << tangentLB.transpose().format(CleanFmt) << " <= z <= " << tangentUB.transpose().format(CleanFmt) << std::endl;

      std::cout << std::endl << "Linear Constraints:" << std::endl;
      std::cout << linCstrLB.format(CleanFmt) << " <= linCstr = " << linCstr.format(CleanFmt) << " <= " << linCstrUB.format(CleanFmt) << std::endl;
      std::cout << "gradlinCstr = " << diffLinCstr.transpose().format(CleanFmt) << std::endl;

      std::cout << std::endl << "NonLinear Constraints:" << std::endl;
      std::cout << nonLinCstrLB.format(CleanFmt) << " <= NonLinCstr = " << nonLinCstr.format(CleanFmt) << " <= " << nonLinCstrUB.format(CleanFmt) << std::endl;
      std::cout << "gradNonLinCstr = " << diffNonLinCstr.transpose().format(CleanFmt) << std::endl;
    }
  };
  struct SolverOption
  {
    int maxIter = 1000;
    double epsilon_P;
    double epsilon_D;
  };
}
#endif //_PGS_TEMPDATA_H_
