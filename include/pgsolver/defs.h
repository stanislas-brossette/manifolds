#ifndef _PGS_DEFS_H_
#define _PGS_DEFS_H_

#include <Eigen/Core>

namespace pgs
{
  typedef Eigen::Ref<Eigen::VectorXd> RefVec;
  typedef Eigen::Ref<const Eigen::VectorXd> ConstRefVec;
  typedef Eigen::Ref<Eigen::MatrixXd> RefMat;
  typedef Eigen::Ref<const Eigen::MatrixXd> ConstRefMat;
  typedef Eigen::Ref<Eigen::VectorXd>::SegmentReturnType Segment;
  typedef Eigen::Ref<const Eigen::VectorXd>::ConstSegmentReturnType ConstSegment;
  typedef Eigen::VectorXd::Index Index;
}

#endif //_PGS_DEFS_H_

