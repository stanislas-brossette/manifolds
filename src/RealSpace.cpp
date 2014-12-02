#include <pgsolver/RealSpace.h>

namespace pgs
{
  RealSpace::RealSpace(Index n)
    : Manifold(n, n)
  {
  }
  
  bool RealSpace::isValidInit(const Eigen::VectorXd& val ) const
  {
    bool out( dim() == val.size());
    return out;
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

  ConstSegment RealSpace::getValueConst(ConstRefVec& val, size_t i) const
  {
    assert(i < 1 && "invalid index");
    return val.segment(0,static_cast<long> (representationDim()));
  }

  std::string RealSpace::toString(ConstRefVec& val, std::string& prefix) const
  {
    std::stringstream ss;
    ss << prefix << val.transpose() << std::endl;
    return ss.str();
  }

  Segment RealSpace::getValue(RefVec val, size_t i) const
  {
    assert(i < 1 && "invalid index");
    return val.segment(0, static_cast<long> (representationDim()));
  }

  void RealSpace::plus_(RefVec out, ConstRefVec& x, ConstRefVec& v) const
  {
    out = x + v;
  }

  void RealSpace::minus_(RefVec out, ConstRefVec& x, ConstRefVec& y) const
  {
    out = x - y;
  }

  void RealSpace::setIdentity_(RefVec out) const
  {
    out.setZero();
  }
}
