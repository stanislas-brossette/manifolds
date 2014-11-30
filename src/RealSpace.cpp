#include <pgsolver/RealSpace.h>

namespace pgs
{
  RealSpace::RealSpace(Index n)
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
    return val.segment(0,static_cast<long> (representationDim()));
  }

  std::string RealSpace::toString(const Eigen::Ref<const Eigen::VectorXd>& val, std::string& prefix) const
  {
    std::stringstream ss;
    ss << prefix << val.transpose() << std::endl;
    return ss.str();
  }

  Segment RealSpace::getValue(Eigen::Ref<Eigen::VectorXd> val, size_t i) const
  {
    assert(i < 1 && "invalid index");
    return val.segment(0, static_cast<long> (representationDim()));
  }

  void RealSpace::plus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v) const
  {
    out = x + v;
  }

  void RealSpace::minus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y) const
  {
    out = x - y;
  }

  void RealSpace::setIdentity_(Eigen::Ref<Eigen::VectorXd> out) const
  {
    out.setZero();
  }
}
