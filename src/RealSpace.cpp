#include "RealSpace.h"

namespace pgs
{
  RealSpace::RealSpace(size_t n)
    : Manifold(n, n)
  {
  }

  size_t RealSpace::numberOfSubmanifolds() const
  {
    return 1;
  }
  const Manifold& RealSpace::operator()(size_t i) const
  {
    assert(i < 1 && "invalid index");
    return *this;
  }

  ConstSegment RealSpace::getValueConst(const Eigen::Ref<const Eigen::VectorXd>& val, size_t i) const
  {
    assert(i < 1 && "invalid index");
    return val.segment(0,representationDim());
  }

  Segment RealSpace::getValue(Eigen::Ref<Eigen::VectorXd> val, size_t i) const
  {
    assert(i < 1 && "invalid index");
    return val.segment(0, representationDim());
  }

  void RealSpace::plus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v) const
  {
    out = x + v;
  }

  void RealSpace::setIdentity_(Eigen::Ref<Eigen::VectorXd> out) const
  {
    out.setZero();
  }
}