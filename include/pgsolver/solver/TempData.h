#ifndef _PGS_TEMPDATA_H_
#define _PGS_TEMPDATA_H_

namespace pgs
{
  struct ProblemEvaluation
  {
    double obj;                           // Value of the costFunction
    Eigen::VectorXd diffObj;              // Value ot the differential of the costFunction  
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
  };
  struct SolverOption
  {
    int maxIter = 1000;
    double epsilon_P;
    double epsilon_D;
  };
}
#endif //_PGS_TEMPDATA_H_
