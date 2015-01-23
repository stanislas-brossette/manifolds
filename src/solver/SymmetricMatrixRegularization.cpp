#include "SymmetricMatrixRegularization.h"

namespace
{
  //Eigen decomposition for a 2d matrix
  void eig2d(const Eigen::Ref<Eigen::Matrix2d>& M, Eigen::Matrix2d& Q, Eigen::Vector2d& d)
  {
    double b = M.trace();           //TODO: ensure numerical robustness to catastrophic cancelation
    double c = M.determinant();
    double delta = b *b - 4 * c;
    assert(delta >= 0);

    d(0) = (b + sqrt(delta)) / 2;
    d(1) = b - d(0);
    Q.coeffRef(0, 0) = M.coeffRef(0, 1);
    Q.coeffRef(1, 0) = d(0) - M.coeffRef(0, 0);
    Q.col(0) /= Q.col(0).norm();
    Q.coeffRef(0, 1) = -Q.coeffRef(1, 0);
    Q.coeffRef(1, 1) = Q.coeffRef(0, 0);
  }
}

namespace pgs
{
  SymmetricMatrixRegularization::SymmetricMatrixRegularization(int size, double valMin)
    :valMin_(valMin), tmp1_(size, size), tmp2_(size,size), decomposition_(size)
  {
  }

  const Eigen::MatrixXd& SymmetricMatrixRegularization::regularize(const Eigen::MatrixXd& M)
  {
    tmp1_ = M;
    regularizeInPlace(tmp1_);
    return tmp1_;
  }

  void SymmetricMatrixRegularization::regularizeInPlace(Eigen::MatrixXd& M)
  {
    assert(M.rows() == M.cols());
    assert(M.rows() == decomposition_.rows());
    decomposition_.compute(M);
    MatrixB B = decomposition_.matrixB();
    MatrixL L = decomposition_.matrixL();

    tmp2_.setIdentity();
    //P
    tmp2_ = decomposition_.transpositionsP()*tmp2_;

    //L'*P
    M.noalias() = L.transpose()*tmp2_;
    
    //Q D Q' L' P
    for (Index i = 0; i < B.blocks(); ++i)
    {
      Index k = B.blockStart(i);
      if (B.blockSize(i) == 1)
      {
        double d = std::max(B.block(i).value(), valMin_);
        tmp2_.row(k) = d*M.row(k);
      }
      else
      {
        Eigen::Matrix2d Q;
        Eigen::Vector2d d;
        eig2d(B.block(i), Q, d);
        d(0) = std::max(d(0), valMin_);
        d(1) = std::max(d(1), valMin_);
        Eigen::Matrix2d DQt = d.asDiagonal()*Q.transpose();
        Eigen::Matrix2d QDQt = Q*DQt;
        tmp2_.middleRows(k, 2).noalias() = QDQt * M.middleRows(k, 2);
      }
    }

    //L Q D Q' L' P
    M = L*tmp2_;

    //P L Q D Q' L' P
    M = decomposition_.transpositionsP()*tmp2_;
  }
}
