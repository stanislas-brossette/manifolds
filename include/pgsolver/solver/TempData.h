#ifndef _PGS_TEMPDATA_H_
#define _PGS_TEMPDATA_H_

#include <Eigen/Core>

namespace pgs
{
  /// \brief Structure containing the current evaluation of the problem.
  /// Its integrity is NOT guarantied. In the sence that all the values are not
  /// necessarily evaluated at the same point.
  struct ProblemEvaluation
  {
    /// \brief Dimension of the variable space (Tangent space of
    /// problems'Manifold
    Index varDim;
    /// \brief total dimension of Linear cstr
    Index linCstrDim;
    /// \brief total dimension of Non-Linear cstr
    Index nonLinCstrDim;


    /// \brief Value of the costFunction
    double obj;
    /// \brief Value ot the differential of the costFunction
    Eigen::MatrixXd diffObj;
    /// \brief Value ot the differential of the costFunction at the previous
    /// iteration
    Eigen::MatrixXd prevDiffObj;
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
    /// \brief Value of the differential of the NonLinear constraints at the
    /// previous iteration
    Eigen::MatrixXd prevDiffNonLinCstr;
    /// \brief Lower Bounds on the NonLinear Constraints
    Eigen::VectorXd nonLinCstrLB;
    /// \brief Upper Bounds on the NonLinear Constraints
    Eigen::VectorXd nonLinCstrUB;

    /// \brief Hessian of the Lagrangian
    Eigen::MatrixXd Hessian;
    /// \brief Hessian of the cost Function
    Eigen::MatrixXd HessianCost;
    /// \brief Vector containing all the Hessians of the non-Linear constraints
    std::vector<Eigen::MatrixXd> HessiansCstr;

    /// \brief Lagrangian
    double lag;
    /// \brief Differential of the Lagrangian
    Eigen::MatrixXd diffLag;
    /// \brief Differential of the Lagrangian at the previous iteration
    Eigen::MatrixXd prevDiffLag;

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

    /// \brief Concatenation of -infLinCstr and -infNonLinCstr
    /// allInfCstr(i) = LB - c(i)
    Eigen::VectorXd allInfCstr;
    /// \brief Concatenation of -supLinCstr and -supNonLinCstr
    /// allSupCstr(i) = UB - c(i)
    Eigen::VectorXd allSupCstr;

    /// \brief Concatenation of linear and nonLinear constraints LB
    Eigen::VectorXd allCstrLB;
    /// \brief Concatenation of linear and nonLinear constraints UB
    Eigen::VectorXd allCstrUB;
    /// \brief Concatenation of linear and nonLinear constraints
    Eigen::VectorXd allCstr;
    /// \brief Concatenation of linear and nonLinear diffConstraints
    Eigen::MatrixXd allDiffCstr;

    //-------------- FEASIBILITY DATA ---------------------
    /// \brief Vector used in the calculation of the cost function of the
    /// feasibility problem
    Eigen::VectorXd feasibilityCostF;
    /// \brief Matrix used for describing linearization of feasibility
    /// constraints in the following order:
    /// [ diffLinConstraint 0     ][1               ][-1              ]
    /// [ diffLinConstraint 1     ][  1             ][  -1            ]
    /// [ diffLinConstraint .     ][    1           ][    -1          ]
    /// [ diffLinConstraint N     ][      1         ][      -1        ]
    /// [ diffNonLinConstraint 0  ][        1       ][        -1      ]
    /// [ diffNonLinConstraint 1  ][          1     ][          -1    ]
    /// [ diffNonLinConstraint .  ][            1   ][            -1  ]
    /// [ diffNonLinConstraint N  ][              1 ][              -1]
    Eigen::MatrixXd feasibilityAllDiffCstr;
    Eigen::VectorXd feasibilityLB;
    Eigen::VectorXd feasibilityUB;
    Eigen::VectorXd feasibleValue;
    Eigen::VectorXd infeasibilityInf;
    Eigen::VectorXd infeasibilitySup;
    Index nFeasCstr;

    /// \brief Vector containing the status of each cstr in the same order as
    /// in allCstr with values taken in enum eCstrStatus
    /// REMINDER:
    /// eCstrStatus {VIOLATED_LB = -2,VIOLATED_UB = -1,SATISFIED = 0, 
    /// ACTIVE_LB = 1, ACTIVE_UB = 2, ACTIVE_EQUALITY = 3}
    Eigen::VectorXi infeasStatus;

    //-------------- RESTORATION DATA ---------------------
    Eigen::VectorXd restorationDiffObj;
    Eigen::MatrixXd restorationAllDiffCstr;
    Eigen::VectorXd restorationAllInfCstr;
    Eigen::VectorXd restorationAllSupCstr;
    Eigen::VectorXd restorationLagMult;


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
}
#endif //_PGS_TEMPDATA_H_
