#ifndef _MANIFOLDS_VIEW_H_
#define _MANIFOLDS_VIEW_H_

#include <Eigen/Core>
#include <manifolds/defs.h>

namespace pgs
{
  enum eDimension
  {
    R,  //representation space
    T,  //tangent space
    F,  //full space
  };

  template<int Dr, int Dc> struct ViewReturnType { typedef Eigen::Block<RefMat> Type; };
  template<int Dr, int Dc> struct ConstViewReturnType { typedef const Eigen::Block<ConstRefMat> Type; };
  template<int Dc> struct ViewReturnType<F, Dc> { typedef RefMat::ColsBlockXpr Type; };
  template<int Dc> struct ConstViewReturnType<F, Dc> { typedef ConstRefMat::ConstColsBlockXpr Type; };
  template<int Dr> struct ViewReturnType<Dr, F> { typedef RefMat::RowsBlockXpr Type; };
  template<int Dr> struct ConstViewReturnType<Dr, F> { typedef ConstRefMat::ConstRowsBlockXpr Type; };
  template<> struct ViewReturnType<F, F> { typedef RefMat Type; };
  template<> struct ConstViewReturnType<F, F> { typedef ConstRefMat Type; };
}

#endif //_MANIFOLDS_VIEW_H_

