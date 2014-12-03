#ifndef _PGS_EXPMAPMATRIX_H_
#define _PGS_EXPMAPMATRIX_H_

#include <Eigen/Core>
#include <pgsolver/defs.h>

namespace pgs
{
  struct ExpMapMatrix
  {
    static const double prec;
    static const int OutputDim_ = 9;
    static const int InputDim_ = 3;
    typedef Eigen::Matrix3d DisplayType;
    typedef Eigen::Matrix<double, 9, 1> OutputType;
    static bool isValidInit(const Eigen::VectorXd& val);
    static void plus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v);
    static void minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& y);
    static void setIdentity_(RefVec out);
    static Eigen::MatrixXd diffMap_(const ConstRefVec& x);
    static void applyDiffMap_(
        RefMat out, const ConstRefMat& in, const ConstRefVec& x);
  };
}

#endif //_PGS_EXPMAPMATRIX_H_
