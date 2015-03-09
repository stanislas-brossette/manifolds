#include <limits>

#include <manifolds/S2.h>
#include <manifolds/pgs_assert.h>

namespace pgs
{
  S2::S2()
    : Manifold(2, 3, 3)
  {
    name() = "S2";
    setTypicalMagnitude(Eigen::Vector3d::Constant(M_PI));
  }

  S2::S2(double magnitude)
    : Manifold(2, 3, 3)
  {
    name() = "S2";
    setTypicalMagnitude(Eigen::Vector3d::Constant( magnitude));
  }
  S2::S2(const ConstRefVec& magnitude)
    : Manifold(2, 3, 3)
  {
    pgs_assert(magnitude.size() == 3 && "magnitude on R^n must be of size n");
    name() = "S2";
    setTypicalMagnitude (magnitude);
  }

  bool S2::isInM_(const Eigen::VectorXd& val ) const
  {
    bool out(val.lpNorm<2>() == 1.0);
    return out;
  }

  size_t S2::numberOfSubmanifolds() const
  {
    return 1;
  }
  const Manifold& S2::operator()(size_t i) const
  {
    pgs_assert(i < 1 && "invalid index");
    return *this;
  }

  std::string S2::toString(const ConstRefVec& val, const std::string& prefix) const
  {
    std::string matPrefix = prefix + '[';
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", matPrefix, "]");
    std::stringstream ss;
    ss << val.transpose().format(CleanFmt);
    return ss.str();
  }

  void S2::retractation_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    Eigen::Vector3d sum;
    sum = x+v;
    out = sum/sum.lpNorm<2>();
  }

  void S2::pseudoLog_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    logarithm(out, x, y);
  }

  void S2::logarithm (RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    Eigen::Vector3d projDiff;
    proj(projDiff,x,y-x);
    out = distance(x,y)*projDiff/projDiff.lpNorm<2>();
  }

  void S2::proj (RefMat out, const ConstRefVec& x, const ConstRefMat& y) const
  {
    out = y - (x.transpose()*y).trace()*x;
  }
    
  double S2::distance (const ConstRefVec& x, const ConstRefVec& y) const
  {
    return acos(x.dot(y));
  }

  void S2::pseudoLog0_(RefVec , const ConstRefVec& ) const
  {
    pgs_assert(1 && "Unimplemented method in S2");
  }

  void S2::setZero_(RefVec ) const
  {
    pgs_assert(1 && "Unimplemented method in S2");
  }

  Eigen::MatrixXd S2::diffRetractation_(const ConstRefVec& out) const
  {
    pgs_assert(1 && "Unimplemented method in S2");
    return out;
  }

  void S2::applyDiffRetractation_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    proj(out, in, x);
  }

  Eigen::MatrixXd S2::diffPseudoLog0_(const ConstRefVec& out) const
  {
    pgs_assert(1 && "Unimplemented method in S2");
    return out;
  }

  void S2::applyDiffPseudoLog0_(RefMat , const ConstRefMat& , const ConstRefVec& ) const
  {
    pgs_assert(1 && "Unimplemented method in S2");
  }

  void S2::applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    Eigen::Vector3d y;
    retractation(y, x, v);
    proj(out, in, y);
  }

  void S2::applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& ) const
  {
    proj(out, in, x);
  }

  void S2::tangentConstraint_(RefMat, const ConstRefVec&) const
  {
  }

  bool S2::isInTxM_(const ConstRefVec& x, const ConstRefVec& v) const
  {
    bool out = (x.dot(v) == 0.0);
    return out;
  }

  void S2::forceOnTxM_(RefVec out, const ConstRefVec& v, const ConstRefVec& x) const
  {
    out = v - x.dot(v)*v;
  }

  void S2::limitMap_(RefVec out) const
  {
    out.setConstant(std::numeric_limits<double>::infinity());
  }

  void S2::getTypicalMagnitude_(RefVec out) const
  {
    out = typicalMagnitude_;
  }

  void S2::setTypicalMagnitude(double magnitude)
  {
    setTypicalMagnitude (Eigen::VectorXd::Constant(tangentDim(), magnitude));
  }

  void S2::setTypicalMagnitude(const ConstRefVec& out)
  {
    typicalMagnitude_ = out;
  }
}

