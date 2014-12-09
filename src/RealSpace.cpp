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

  //Segment RealSpace::getValue(RefVec val, size_t i) const
  //{
  //  assert(val.size() == representationDim());
  //  assert(i < 1 && "invalid index");
  //  return val.segment(0, static_cast<long> (representationDim()));
  //}

  //ConstSegment RealSpace::getValueConst(const ConstRefVec& val, size_t i) const
  //{
  //  assert(val.size() == representationDim());
  //  assert(i < 1 && "invalid index");
  //  return val.segment(0,static_cast<long> (representationDim()));
  //}

  //Segment RealSpace::getValueTangent(RefVec val, size_t i) const
  //{
  //  assert(val.size() == dim());
  //  assert(i < 1 && "invalid index");
  //  return val.segment(0, static_cast<long> (dim()));
  //}

  //ConstSegment RealSpace::getValueTangentConst(const ConstRefVec& val, size_t i) const
  //{
  //  assert(val.size() == dim());
  //  assert(i < 1 && "invalid index");
  //  return val.segment(0,static_cast<long> (dim()));
  //}

  std::string RealSpace::toString(const ConstRefVec& val, const std::string& prefix) const
  {
    std::string matPrefix = prefix + '[';
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", matPrefix, "]");
    std::stringstream ss;
    ss << val.transpose().format(CleanFmt);
    return ss.str();
  }

  void RealSpace::plus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    out = x + v;
  }

  void RealSpace::minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    out = x - y;
  }

  void RealSpace::invMap_(RefVec out, const ConstRefVec& x) const
  {
    out = x;
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
