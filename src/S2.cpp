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

  bool S2::isInM_(const Eigen::VectorXd& val , const double& prec) const
  {
    bool out(fabs(val.lpNorm<2>() - 1.0) < prec);
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

  void S2::createRandomPoint_(RefVec out, double ) const
  {
    rand(out);
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
    Eigen::Vector3d diff, projDiff;
    diff = y-x;
    projVec(projDiff, diff, x);
    out = distance(x,y)*projDiff/projDiff.lpNorm<2>();
  }

  double S2::distance (const ConstRefVec& x, const ConstRefVec& y) const
  {
    return acos(x.dot(y));
  }

  void S2::projRows(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    for (Index i = 0; i < in.rows(); ++i)
      out.row(i) = in.row(i) - (x*in.row(i)).trace()*x.transpose();
  }
  void S2::projCols(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    for (Index i = 0; i < in.cols(); ++i)
      projVec(out.col(i), in.col(i), x);
  }
  void S2::projVec(RefVec out, const ConstRefVec& in, const ConstRefVec& x) const
  {
      out = in - (x.dot(in))*x;
  }

  void S2::rand(RefVec out) const
  {
    pgs_assert( out.size() == 3 && "Wrong size of output for S2::rand");
    out = Eigen::Vector3d::Random();
    double nout = out.lpNorm<2>();
    out = out/nout;
  }

  void S2::randVec(RefVec out, const ConstRefVec& x) const
  {
    pgs_assert( out.size() == 3 && "Wrong size of output for S2::rand");
    pgs_assert( x.size() == 3 && "Wrong size of output for S2::rand");
    out = Eigen::Vector3d::Random();
    projVec(out, out, x);
    double nout = out.lpNorm<2>();
    out = out/nout;
  }
  Eigen::Vector3d S2::randVec(const ConstRefVec& x) const
  {
    Eigen::Vector3d out;
    randVec(out, x);
    return out;
  }
    
  void S2::pseudoLog0_(RefVec , const ConstRefVec& ) const
  {
    throw std::runtime_error("Unimplemented method in S2");
  }

  void S2::setZero_(RefVec ) const
  {
    throw std::runtime_error("Unimplemented method in S2");
  }

  Eigen::MatrixXd S2::diffRetractation_(const ConstRefVec& out) const
  {
    throw std::runtime_error("Unimplemented method in S2");
    return out;
  }

  void S2::applyDiffRetractation_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    projRows(out, in, x);
  }

  Eigen::MatrixXd S2::diffPseudoLog0_(const ConstRefVec& out) const
  {
    throw std::runtime_error("Unimplemented method in S2");
    return out;
  }

  void S2::applyDiffPseudoLog0_(RefMat , const ConstRefMat& , const ConstRefVec& ) const
  {
    throw std::runtime_error("Unimplemented method in S2");
  }

  void S2::applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    Eigen::Vector3d y;
    retractation(y, x, v);
    projCols(out, in, y);
  }

  void S2::applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& ) const
  {
    projRows(out, in, x);
  }

  void S2::tangentConstraint_(RefMat out, const ConstRefVec& x) const
  {
    out = x.transpose();
  }

  bool S2::isInTxM_(const ConstRefVec& x, const ConstRefVec& v, const double& prec) const
  {
    bool out = (fabs(x.dot(v)) < prec);
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

