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
    static void plus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v);
    static void minus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y);
    static void setIdentity_(Eigen::Ref<Eigen::VectorXd> out);
  };
}

#endif //_PGS_EXPMAPMATRIX_H_
