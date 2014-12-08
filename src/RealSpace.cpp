#include <pgsolver/RealSpace.h>

namespace pgs
{
  RealSpace::RealSpace(Index n)
    : Manifold(n, n)
  {
  }
  
  bool RealSpace::isValidInit_(const Eigen::VectorXd& val ) const
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

  ConstSegment RealSpace::getValueConst(const ConstRefVec& val, size_t i) const
  {
    assert(i < 1 && "invalid index");
    return val.segment(0,static_cast<long> (representationDim()));
  }

  std::string RealSpace::toString(const ConstRefVec& val, std::string& prefix) const
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

  void RealSpace::plus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    out = x + v;
  }

  void RealSpace::minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    out = x - y;
  }

  void RealSpace::setIdentity_(RefVec out) const
  {
    out.setZero();
  }

  Eigen::MatrixXd RealSpace::diffMap_(const ConstRefVec& ) const
  {
    return Eigen::MatrixXd::Identity(representationDim(),dim());
  }

  void RealSpace::applyDiffMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& ) const
  {
    out = in;
  }
  
  Eigen::MatrixXd RealSpace::diffInvMap_(const ConstRefVec&) const
  {
    return Eigen::MatrixXd::Identity(representationDim(),dim());
  }

  void RealSpace::applyDiffInvMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& ) const
  {
    out = in;
  }
}
