#ifndef _PGS_DEFS_H_
#define _PGS_DEFS_H_

#include <Eigen/Core>

namespace pgs
{
  typedef Eigen::Ref<Eigen::VectorXd>::SegmentReturnType Segment;
  typedef Eigen::Ref<const Eigen::VectorXd>::ConstSegmentReturnType ConstSegment;
  typedef Eigen::VectorXd::Index Index;
}

#endif //_PGS_DEFS_H_

