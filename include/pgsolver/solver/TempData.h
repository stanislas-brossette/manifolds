#ifndef _PGS_TEMPDATA_H_
#define _PGS_TEMPDATA_H_

#include <Eigen/Core>

namespace pgs
{
  struct ProblemEvaluation
  {
    /// \brief Value of the costFunction
    double obj;                           
    /// \brief Value ot the differential of the costFunction  
    Eigen::MatrixXd diffObj;              
    /// \brief Tangent Lower Bound
    Eigen::VectorXd tangentLB;            
    /// \brief Tangent Upper Bound
    Eigen::VectorXd tangentUB;            

    /// \brief Value of the Linear Constraints
    Eigen::VectorXd linCstr;              
    /// \brief Value of the differential of the Linear constraints
    Eigen::MatrixXd diffLinCstr;          
    /// \brief Lower Bounds on the Linear Constraints 
    Eigen::VectorXd linCstrLB;            
    /// \brief Upper Bounds on the Linear Constraints
    Eigen::VectorXd linCstrUB;            

    /// \brief Value of the NonLinear Constraints
    Eigen::VectorXd nonLinCstr;           
    /// \brief Value of the differential of the NonLinear constraints
    Eigen::MatrixXd diffNonLinCstr;       
    /// \brief Lower Bounds on the NonLinear Constraints 
    Eigen::VectorXd nonLinCstrLB;         
    /// \brief Upper Bounds on the NonLinear Constraints
    Eigen::VectorXd nonLinCstrUB;         
    
    /// \brief Hessian of the Lagrangian
    Eigen::MatrixXd Hessian;              

    /// \brief Lagrangian
    double lag;                           
    /// \brief Differential of the Lagrangian
    Eigen::MatrixXd diffLag;              

    /// \brief linearized inf bound of linear constraints
    /// LinConstraint - infBound
    Eigen::VectorXd linearizedInfBndLinCstr;
    /// \brief linearized inf bound of nonlinear constraints 
    /// nonLinConstraint - infBound
    Eigen::VectorXd linearizedInfBndNonLinCstr;
    /// \brief linearized sup bound of linear constraints (supBound -
    /// LinConstraint - supBound
    Eigen::VectorXd linearizedSupBndLinCstr;
    /// \brief linearized sup bound of nonlinear constraints (supBound -
    /// nonLinConstraint - supBound
    Eigen::VectorXd linearizedSupBndNonLinCstr;

    void print()
    {
      Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
      std::cout << std::endl << "Problem evaluation" << std::endl;

      std::cout << std::endl << "Objective Function:" << std::endl;
      std::cout << "f(phi_x(z))=" << obj << std::endl;
      std::cout << "grad_z(f(phi_x(z)):" << std::endl << diffObj.transpose().format(CleanFmt) << std::endl;

      std::cout << std::endl << "Bounds:" << std::endl;
      std::cout << tangentLB.transpose().format(CleanFmt) << " <= z <= " << tangentUB.transpose().format(CleanFmt) << std::endl;

      std::cout << std::endl << "Linear Constraints:" << std::endl;
      for(Index i = 0; i<linCstr.size(); ++i)
      {
        std::cout << "linCstr" << i << " : " << linCstrLB[i] << " <= " << linCstr[i] << " <= " << linCstrUB[i] << std::endl;
      }
      std::cout << "gradlinCstr: " << std::endl << diffLinCstr.transpose().format(CleanFmt) << std::endl;

      std::cout << std::endl << "NonLinear Constraints:" << std::endl;
      for(Index i = 0; i<nonLinCstr.size(); ++i)
      {
        std::cout << "nonLinCstr" << i << " : " << nonLinCstrLB[i] << " <= " << nonLinCstr[i] << " <= " << nonLinCstrUB[i] << std::endl;
      }
      std::cout << "gradNonLinCstr: " << std::endl << diffNonLinCstr.transpose().format(CleanFmt) << std::endl;

      std::cout << std::endl << "Lagrangian:" << std::endl;
      std::cout << "Lag(phi_x(z))=" << lag << std::endl;
      std::cout << "grad_z(Lag(phi_x(z)):" << std::endl << diffLag.transpose().format(CleanFmt) << std::endl;
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
