#ifndef _PGS_SYMETRIC_MATRIX_REGULARIZATION_
#define _PGS_SYMETRIC_MATRIX_REGULARIZATION_

#include <Eigen/Eigen>
#include "BlockDiagonalMatrix.h"
#include "LBLT.h"

namespace pgs
{
  class SymmetricMatrixRegularization
  {
  public:
    typedef Eigen::internal::LBLT_Traits<Eigen::MatrixXd, Eigen::Lower>::MatrixB MatrixB;
    typedef Eigen::internal::LBLT_Traits<Eigen::MatrixXd, Eigen::Lower>::MatrixL MatrixL;
    typedef Eigen::MatrixXd::Index Index;

    SymmetricMatrixRegularization(int size, double valMin);

    const Eigen::MatrixXd& regularize(const Eigen::MatrixXd& M);
    void regularizeInPlace(Eigen::MatrixXd& M);

  private:
    double valMin_;
    Eigen::MatrixXd tmp1_;
    Eigen::MatrixXd tmp2_;
    Eigen::LBLT<Eigen::MatrixXd, Eigen::Lower> decomposition_;
  };
}

#endif //_PGS_SYMETRIC_MATRIX_REGULARIZATION_