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
    /// \brief Differential of the Lagrangian at the previous iteration
    Eigen::MatrixXd prevDiffLag;

    ///// \brief linearized inf bound of bound constraints
    ///// z - infBound
    //Eigen::VectorXd infBndCstr;
    ///// \brief linearized sup bound of bound constraints
    ///// z - supBound
    //Eigen::VectorXd supBndCstr;
    /// \brief linearized inf bound of linear constraints
    /// LinConstraint - infBound
    Eigen::VectorXd infLinCstr;
    /// \brief linearized inf bound of nonlinear constraints 
    /// nonLinConstraint - infBound
    Eigen::VectorXd infNonLinCstr;
    /// \brief linearized sup bound of linear constraints (supBound -
    /// LinConstraint - supBound
    Eigen::VectorXd supLinCstr;
    /// \brief linearized sup bound of nonlinear constraints (supBound -
    /// nonLinConstraint - supBound
    Eigen::VectorXd supNonLinCstr;

    /// \brief Concatenation of infLinCstr and infNonLinCstr 
    Eigen::VectorXd allInfCstr;
    /// \brief Concatenation of supLinCstr and supNonLinCstr 
    Eigen::VectorXd allSupCstr;

    /// \brief Concatenation of linear and nonLinear constraints
    Eigen::VectorXd allCstr; 
    /// \brief Concatenation of linear and nonLinear diffConstraints
    Eigen::MatrixXd allDiffCstr; 

    void print()
    {
      Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");

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

      //std::cout << "AllCstr:" << std::endl << allCstr << std::endl;
      std::cout << std::endl << "Lagrangian:" << std::endl;
      std::cout << "Lag(phi_x(z))=" << lag << std::endl;
      std::cout << "grad_z(Lag(phi_x(z)):" << std::endl << diffLag.transpose().format(CleanFmt) << std::endl;
      std::cout << "Hessian: " << std::endl << Hessian.format(CleanFmt) << std::endl;
    }
  };

  struct LagrangeMultipliers
  {
    Eigen::VectorXd bounds;
    Eigen::VectorXd linear;
    Eigen::VectorXd nonLinear;
    Eigen::VectorXd all;

    void initOnes()
    {
      bounds.setOnes();
      linear.setOnes();
      nonLinear.setOnes();
      all.setZero();
      update(bounds, linear, nonLinear);
    }

    void update(const Eigen::VectorXd& bnd,
                const Eigen::VectorXd& lin, 
                const Eigen::VectorXd& nonLin)
    {
      assert(bnd.size() == bounds.size());
      assert(lin.size() == linear.size());
      assert(nonLin.size() == nonLinear.size());

      bounds = bnd;
      linear = lin;
      nonLinear = nonLin;
      all.head(bounds.size()) = bnd;
      all.segment(bounds.size(), linear.size()) = lin;
      all.tail(nonLinear.size()) = nonLin;
    }
    void print()
    {
      Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
      std::cout << "Lagrange Multipliers Bounds Cstr" << std::endl;
      std::cout << bounds.format(CleanFmt)<< std::endl;
      std::cout << "Lagrange Multipliers Linear Cstr" << std::endl;
      std::cout << linear.format(CleanFmt)<< std::endl;
      std::cout << "Lagrange Multipliers NonLinear Cstr" << std::endl;
      std::cout << nonLinear.format(CleanFmt)<< std::endl;
      //std::cout << "Lagrange Multipliers All Cstr" << std::endl;
      //std::cout << all.format(CleanFmt)<< std::endl;
    }
  };

  struct SolverOption
  {
    int maxIter = 1000;
    double epsilon_P = 1e-6;
    double epsilon_D = 1e-6;
  };

}
#endif //_PGS_TEMPDATA_H_
