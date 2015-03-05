#include <manifolds/RealSpace.h>
#include <manifolds/pgs_assert.h>

namespace pgs
{
  RealSpace::RealSpace(Index n)
    : Manifold(n, n, n)
  {
  }

  bool RealSpace::isInM_(const Eigen::VectorXd& val ) const
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
    pgs_assert(i < 1 && "invalid index");
    return *this;
  }

  std::string RealSpace::toString(const ConstRefVec& val, const std::string& prefix) const
  {
    std::string matPrefix = prefix + '[';
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", matPrefix, "]");
    std::stringstream ss;
    ss << val.transpose().format(CleanFmt);
    return ss.str();
  }

  void RealSpace::retractation_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    out = x + v;
  }

  void RealSpace::pseudoLog_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    out = x - y;
  }

  void RealSpace::invMap_(RefVec out, const ConstRefVec& x) const
  {
    out = x;
  }

  void RealSpace::setZero_(RefVec out) const
  {
    out.setZero();
  }

  Eigen::MatrixXd RealSpace::diffRetractation_(const ConstRefVec& ) const
  {
    return Eigen::MatrixXd::Identity(representationDim(),dim());
  }

  void RealSpace::applyDiffRetractation_(RefMat out, const ConstRefMat& in, const ConstRefVec& ) const
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

  void RealSpace::applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec& ) const
  {
    out = in;
  }

  void RealSpace::applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec& ) const
  {
    out = in;
  }

  void RealSpace::tangentConstraint_(RefMat, const ConstRefVec&) const
  {
    //matrix is 0xt, no need to fill it.
  }

  bool RealSpace::isInTxM_(const ConstRefVec&, const ConstRefVec&) const
  {
    return true;
  }

  void RealSpace::forceOnTxM_(RefVec out, const ConstRefVec& in, const ConstRefVec&) const
  {
    out = in;
  }
}
