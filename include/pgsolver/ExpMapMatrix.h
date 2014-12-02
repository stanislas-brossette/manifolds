#ifndef _PGS_EXPMAPMATRIX_H_
#define _PGS_EXPMAPMATRIX_H_

#include <Eigen/Core>

namespace pgs
{
  struct ExpMapMatrix
  {
    static int OutputDim();
    typedef Eigen::Matrix3d DisplayType;
    typedef Eigen::VectorXd OutputType;
    static void plus_(RefVec out, ConstRefVec& x, ConstRefVec& v);
    static void minus_(RefVec out, ConstRefVec& x, ConstRefVec& y);
    static void setIdentity_(RefVec out);
    static bool isValidInit(const Eigen::VectorXd& val);
  };
}

#endif //_PGS_EXPMAPMATRIX_H_
