#include <limits>

#include <manifolds/RealSpace.h>
#include <manifolds/pgs_assert.h>

namespace pgs
{
  RealSpace::RealSpace(Index n)
    : Manifold(n, n, n)
  {
    name() = "R" + std::to_string( n );
    typicalMagnitude_.resize(n);
    Eigen::VectorXd temp = Eigen::VectorXd::Constant(n, 1.0);
    setTypicalMagnitude(temp);
  }

  RealSpace::RealSpace(Index n, double magnitude)
    : Manifold(n, n, n)
  {
    name() = "R" + std::to_string( n );
    typicalMagnitude_.resize(n);
    Eigen::VectorXd temp = Eigen::VectorXd::Constant(n, magnitude);
    setTypicalMagnitude(temp);
  }
  RealSpace::RealSpace(Index n, const ConstRefVec& magnitude)
    : Manifold(n, n, n)
  {
    pgs_assert(magnitude.size() == n && "magnitude on R^n must be of size n");
    name() = "R" + std::to_string( n );
    typicalMagnitude_.resize(n);
    setTypicalMagnitude (magnitude);
  }

  bool RealSpace::isInM_(const Eigen::VectorXd& val , const double& ) const
  {
    bool out( dim() == val.size());
    return out;
  }

  void RealSpace::forceOnM_(RefVec out, const ConstRefVec& in) const
  {
    out = in;
  }

  size_t RealSpace::numberOfSubmanifolds() const
  {
    return 1;
  }

  bool RealSpace::isElementary() const
  {
    return true;
  }

  const Manifold& RealSpace::operator()(size_t i) const
  {
    pgs_assert(i < 1 && "invalid index");
    return *this;
  }

  std::string RealSpace::toString(const ConstRefVec& val, const std::string& prefix, int prec) const
  {
    std::string matPrefix = prefix + '[';
    Eigen::IOFormat CleanFmt(prec, 0, ", ", "\n", matPrefix, "]");
    std::stringstream ss;
    ss << val.transpose().format(CleanFmt);
    return ss.str();
  }

  void RealSpace::createRandomPoint_(RefVec out, double coeff) const
  {
    out = coeff*Eigen::VectorXd::Random(representationDim());
  }

  void RealSpace::retractation_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    out = x + v;
  }

  void RealSpace::pseudoLog_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    out = y - x;
  }

  void RealSpace::pseudoLog0_(RefVec out, const ConstRefVec& x) const
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

  Eigen::MatrixXd RealSpace::diffPseudoLog0_(const ConstRefVec&) const
  {
    return Eigen::MatrixXd::Identity(representationDim(),dim());
  }

  void RealSpace::applyDiffPseudoLog0_(RefMat out, const ConstRefMat& in, const ConstRefVec& ) const
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

  bool RealSpace::isInTxM_(const ConstRefVec&, const ConstRefVec&, const double& ) const
  {
    return true;
  }

  void RealSpace::forceOnTxM_(RefVec out, const ConstRefVec& in, const ConstRefVec&) const
  {
    out = in;
  }

  void RealSpace::limitMap_(RefVec out) const
  {
    out.setConstant(std::numeric_limits<double>::infinity());
  }

  void RealSpace::getTypicalMagnitude_(RefVec out) const
  {
    out = typicalMagnitude_;
  }

  void RealSpace::setTypicalMagnitude(double magnitude)
  {
    Eigen::VectorXd temp = Eigen::VectorXd::Constant(tangentDim(), magnitude);
    setTypicalMagnitude(temp);
  }

  void RealSpace::setTypicalMagnitude(const ConstRefVec& out)
  {
    typicalMagnitude_ = out;
  }

  long RealSpace::getTypeId()
  {
    long typeId = utils::hash::computeHash("RealSpace");
    return typeId;
  }
}
