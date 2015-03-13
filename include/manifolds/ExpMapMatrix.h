#ifndef _MANIFOLDS_EXPMAPMATRIX_H_
#define _MANIFOLDS_EXPMAPMATRIX_H_

#include <Eigen/Core>
#include <manifolds/defs.h>
#include <manifolds/ReusableTemporaryMap.h>

namespace pgs
{
  /// \brief Structure representing the exponential map going from
  /// \f$ \mathbb{R}^3 \f$ to SO(3) represented in \f$ \mathbb{R}^{3\times3} \f$
  struct MANIFOLDS_API ExpMapMatrix
  {
    /// \brief precision constant
    static const double prec;

    /// \brief dimension of \f$ \mathbb{R}^{3\times3}=9 \f$ 
    static const int OutputDim_ = 9;

    /// \brief dimension of \f$ \mathbb{R}^3=3 \f$ 
    static const int InputDim_ = 3;
    typedef Eigen::Matrix3d DisplayType;
    typedef Eigen::Matrix3d OutputType;
    static bool isInM_(const Eigen::VectorXd& val, const double& prec);
    static void forceOnM_(RefVec out, const ConstRefVec& in);
    static void retractation_(RefVec out, const ConstRefVec& x, const ConstRefVec& v);
    static void pseudoLog_(RefVec out, const ConstRefVec& x, const ConstRefVec& y);
    static void pseudoLog0_(RefVec out, const ConstRefVec& x);
    static void setZero_(RefVec out);

    static void logarithm(RefVec out, const OutputType& M);
    static void exponential(OutputType& out, const ConstRefVec& v);

    static Eigen::Matrix<double, 9, 3> diffRetractation_(const ConstRefVec& x);
    static void applyDiffRetractation_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, ReusableTemporaryMap& m);
    static Eigen::Matrix<double, 3, 9> diffPseudoLog0_(const ConstRefVec& x);
    static void applyDiffPseudoLog0_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, ReusableTemporaryMap& m);
    static void applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v, ReusableTemporaryMap& m);
    static void applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v, ReusableTemporaryMap& m);

    static void tangentConstraint_(RefMat out, const ConstRefVec& x);
    static bool isInTxM_(const ConstRefVec& x, const ConstRefVec& v, const double& prec);
    static void forceOnTxM_(RefVec out, const ConstRefVec& in, const ConstRefVec& x);

  };
}

#endif //_MANIFOLDS_EXPMAPMATRIX_H_

