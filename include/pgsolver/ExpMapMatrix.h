#ifndef _PGS_EXPMAPMATRIX_H_
#define _PGS_EXPMAPMATRIX_H_

#include <Eigen/Core>
#include <pgsolver/defs.h>
#include <pgsolver/ReusableTemporaryMap.h>

namespace pgs
{
  struct ExpMapMatrix
  {
    static const double prec;
    static const int OutputDim_ = 9;
    static const int InputDim_ = 3;
    typedef Eigen::Matrix3d DisplayType;
    typedef Eigen::Matrix<double, 9, 1> OutputType;
    static bool isValidInit_(const Eigen::VectorXd& val);
    static void plus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v);
    static void minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& y);
    static void invMap_(RefVec out, const ConstRefVec& x);
    static void setIdentity_(RefVec out);

    static void logarithm(RefVec out, const DisplayType& M);
    static void exponential(DisplayType& out, const ConstRefVec& v);

    static Eigen::Matrix<double, 9, 3> diffMap_(const ConstRefVec& x);
    static void applyDiffMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, ReusableTemporaryMap& m);
    static Eigen::Matrix<double, 3, 9> diffInvMap_(const ConstRefVec& x);
    static void applyDiffInvMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, ReusableTemporaryMap& m);
    static void applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v, ReusableTemporaryMap& m);
    static void applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v, ReusableTemporaryMap& m);
  };
}

#endif //_PGS_EXPMAPMATRIX_H_
