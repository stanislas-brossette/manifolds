// Copyright (c) 2015 CNRS
// Authors: Stanislas Brossette, Adrien Escande

// This file is part of manifolds
// manifolds is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.

// manifolds is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// manifolds. If not, see
// <http://www.gnu.org/licenses/>.

#include <limits>
#define _USE_MATH_DEFINES
#include <math.h>

#include <manifolds/utils.h>
#include <manifolds/S2.h>
#include <manifolds/mnf_assert.h>

namespace mnf
{
S2::S2() : Manifold(2, 3, 3)
{
  std::string name("S2");
  setName(name);
  setTypicalMagnitude(Eigen::Vector3d::Constant(M_PI));
  setTrustMagnitude(Eigen::Vector3d::Constant(M_PI));
}

S2::S2(double magnitude) : Manifold(2, 3, 3)
{
  std::string name("S2");
  setName(name);
  setTypicalMagnitude(Eigen::Vector3d::Constant(magnitude));
  setTrustMagnitude(Eigen::Vector3d::Constant(magnitude));
}
S2::S2(const ConstRefVec& magnitude) : Manifold(2, 3, 3)
{
  mnf_assert(magnitude.size() == 3 && "magnitude on R^n must be of size n");
  std::string name("S2");
  setName(name);
  setTypicalMagnitude(magnitude);
  setTrustMagnitude(magnitude);
}

bool S2::isInM_(const Eigen::VectorXd& val, double prec) const
{
  bool out(fabs(val.norm() - 1.0) < prec);
  return out;
}

void S2::forceOnM_(RefVec out, const ConstRefVec& in) const
{
  out = in / in.norm();
}

void S2::getIdentityOnTxM_(RefMat out, const ConstRefVec& x) const
{
  out.setIdentity();
  Eigen::Matrix3d xxt;
  xxt.noalias() = x * x.transpose();
  out -= xxt;
}

size_t S2::numberOfSubManifolds() const { return 1; }

bool S2::isElementary() const { return true; }

const Manifold& S2::operator()(size_t i) const
{
  mnf_assert(i < 1 && "invalid index");
  return *this;
}

std::string S2::toString(const ConstRefVec& val, const std::string& prefix,
                         const Eigen::IOFormat& fmt) const
{
  std::stringstream ss;
  ss << prefix << val.transpose().format(fmt);
  return ss.str();
}

void S2::createRandomPoint_(RefVec out, double) const { rand(out); }

void S2::retractation_(RefVec out, const ConstRefVec& x,
                       const ConstRefVec& v) const
{
  Eigen::Vector3d sum;
  sum = x + v;
  out = sum / sum.norm();
}

void S2::pseudoLog_(RefVec out, const ConstRefVec& x,
                    const ConstRefVec& y) const
{
  logarithm(out, x, y);
}

void S2::logarithm(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
{
  Eigen::Vector3d diff, projDiff;
  diff = y - x;
  projVec(projDiff, diff, x);
  out = distance_(x, y) * projDiff / projDiff.norm();
}

double S2::distance_(const ConstRefVec& x, const ConstRefVec& y) const
{
  // double d = x.dot(y);
  // double res = acos(d);
  double res = 1 - x.dot(y);
  return res;
}

double S2::squaredDistance_(const ConstRefVec& x, const ConstRefVec& y) const
{
  // double d = x.dot(y);
  // double res = pow(acos(d),2);
  double res = pow(1 - x.dot(y), 2);
  return res;
}

double S2::squaredDistance_(const ConstRefVec& x, const ConstRefVec& y, const ConstRefVec& w) const
{
  // double d = x.dot(y);
  // double res = pow(acos(d),2);
  // WARNING: This is ugly. The weight shoud multiply the log. But weight per direction of the tangent space does not make much sense except for RealSpace. So it does not really matter here.
  double res = pow(w(0)*(1 - x.dot(y)), 2);
  return res;
}

Eigen::MatrixXd S2::derivDistanceX_(const ConstRefVec& /*x*/,
                                    const ConstRefVec& y) const
{
  Eigen::Matrix<double, 1, 3> J;

  // double coeff = -1/sqrt(1-pow(x.dot(y),2));
  // J = (coeff*y).transpose();

  J = -y;

  return J;
}
Eigen::MatrixXd S2::derivDistanceY_(const ConstRefVec& x,
                                    const ConstRefVec& /*y*/) const
{
  Eigen::Matrix<double, 1, 3> J;

  // double coeff = -1/sqrt(1-pow(x.dot(y),2));
  // J = (coeff*x).transpose();

  J = -x;

  return J;
}

Eigen::MatrixXd S2::derivSquaredDistanceX_(const ConstRefVec& x,
                                           const ConstRefVec& y) const
{
  Eigen::Matrix<double, 1, 3> J;

  // double coeff = -2*distance(x,y)/sqrt(1-pow(x.dot(y),2));
  // J = (coeff*y).transpose();

  double coeff = 2 * (1 - x.dot(y));
  J = -coeff * y;

  return J;
}
Eigen::MatrixXd S2::derivSquaredDistanceY_(const ConstRefVec& x,
                                           const ConstRefVec& y) const
{
  Eigen::Matrix<double, 1, 3> J;

  // double coeff = -2*distance(x,y)/sqrt(1-pow(x.dot(y),2));
  // J = (coeff*x).transpose();

  double coeff = 2 * (1 - x.dot(y));
  J = -coeff * x;

  return J;
}

Eigen::MatrixXd S2::derivSquaredDistanceX_(const ConstRefVec& x,
                                           const ConstRefVec& y,
                                           const ConstRefVec& w) const
{
  Eigen::Matrix<double, 1, 3> J;

  // double coeff = -2*distance(x,y)/sqrt(1-pow(x.dot(y),2));
  // J = (coeff*y).transpose();

  double coeff = 2 * (1 - x.dot(y));
  // WARNING: This is ugly. The weight shoud multiply the log. But weight per direction of the tangent space does not make much sense except for RealSpace. So it does not really matter here.
  J = -coeff * y * w(0) * w(0);

  return J;
}
Eigen::MatrixXd S2::derivSquaredDistanceY_(const ConstRefVec& x,
                                           const ConstRefVec& y,
                                           const ConstRefVec& w) const
{
  Eigen::Matrix<double, 1, 3> J;

  // double coeff = -2*distance(x,y)/sqrt(1-pow(x.dot(y),2));
  // J = (coeff*x).transpose();

  double coeff = 2 * (1 - x.dot(y));
  // WARNING: This is ugly. The weight shoud multiply the log. But weight per direction of the tangent space does not make much sense except for RealSpace. So it does not really matter here.
  J = -coeff * x * w(0) * w(0);

  return J;
}

void S2::projRows(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
{
  for (Index i = 0; i < in.rows(); ++i)
    out.row(i) = in.row(i) - (x * in.row(i)).trace() * x.transpose();
}
void S2::projCols(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
{
  for (Index i = 0; i < in.cols(); ++i) projVec(out.col(i), in.col(i), x);
}
void S2::projVec(RefVec out, const ConstRefVec& in, const ConstRefVec& x) const
{
  out = in - (x.dot(in)) * x;
}

void S2::rand(RefVec out) const
{
  mnf_assert(out.size() == 3 && "Wrong size of output for S2::rand");
  out = Eigen::Vector3d::Random();
  double nout = out.norm();
  out = out / nout;
}

void S2::randVec(RefVec out, const ConstRefVec& x) const
{
  mnf_assert(out.size() == 3 && "Wrong size of output for S2::rand");
  mnf_assert(x.size() == 3 && "Wrong size of output for S2::rand");
  out = Eigen::Vector3d::Random();
  projVec(out, out, x);
  double nout = out.norm();
  out = out / nout;
}
Eigen::Vector3d S2::randVec(const ConstRefVec& x) const
{
  Eigen::Vector3d out;
  randVec(out, x);
  return out;
}

void S2::pseudoLog0_(RefVec, const ConstRefVec&) const
{
  throw std::runtime_error(
      "Unimplemented method in S2 pseudoLog0_ makes no sense because there is "
      "no 0");
}

void S2::setZero_(RefVec out) const
{
  std::cout << "Warning, you are using SetZero on S2. Returning an arbitrary "
               "point on S2, which is (1, 0, 0)" << std::endl;
  out << 1.0, 0.0, 0.0;
}

Eigen::MatrixXd S2::diffRetractation_(const ConstRefVec& x) const
{
  Eigen::Matrix3d res;
  res << 1-x(0)*x(0), -x(0)*x(1), -x(0)*x(2),
         -x(0)*x(1), 1-x(1)*x(1), -x(1)*x(2),
         -x(0)*x(2), -x(1)*x(2), 1-x(2)*x(2);
  return res;
}

void S2::applyDiffRetractation_(RefMat out, const ConstRefMat& in,
                                const ConstRefVec& x) const
{
  projRows(out, in, x);
}

Eigen::MatrixXd S2::diffPseudoLog0_(const ConstRefVec& out) const
{
  throw std::runtime_error("Unimplemented method in S2");
  return out;
}

void S2::applyDiffPseudoLog0_(RefMat, const ConstRefMat&,
                              const ConstRefVec&) const
{
  throw std::runtime_error("Unimplemented method in S2");
}

void S2::applyTransport_(RefMat out, const ConstRefMat& in,
                         const ConstRefVec& x, const ConstRefVec& v) const
{
  Eigen::Vector3d y;
  retractation(y, x, v);
  Eigen::Matrix3d R = utils::computeRotBetweenVec(x, y);
  Eigen::Vector3d tmp;
  for (Index i = 0; i < in.cols(); ++i)
  {
    tmp = R * in.col(i);
    out.col(i) = tmp;
  }
}

void S2::applyInvTransport_(RefMat out, const ConstRefMat& in,
                            const ConstRefVec& x, const ConstRefVec& v) const
{
  Eigen::Vector3d y;
  retractation(y, x, v);
  Eigen::Matrix3d R = utils::computeRotBetweenVec(x, y);
  Eigen::Vector3d tmp;
  for (Index i = 0; i < in.cols(); ++i)
  {
    tmp = R.transpose() * in.col(i);
    out.col(i) = tmp;
  }
}

void S2::applyInvTransportOnTheRight_(RefMat out, const ConstRefMat& in,
                                      const ConstRefVec& x,
                                      const ConstRefVec& v) const
{
  Eigen::Vector3d y;
  retractation(y, x, v);
  Eigen::Matrix3d R = utils::computeRotBetweenVec(x, y);
  Eigen::Matrix<double, 3, 1> tmp;
  for (Index i = 0; i < in.rows(); ++i)
  {
    tmp = in.row(i) * R.transpose();
    out.row(i) = tmp;
  }
}

void S2::tangentConstraint_(RefMat out, const ConstRefVec& x) const
{
  out = x.transpose();
}

bool S2::isInTxM_(const ConstRefVec& x, const ConstRefVec& v,
                  const double& prec) const
{
  bool out = (fabs(x.dot(v)) < prec);
  return out;
}

void S2::forceOnTxM_(RefVec out, const ConstRefVec& v,
                     const ConstRefVec& x) const
{
  out = v - x.dot(v) * x;
}

void S2::limitMap_(RefVec out) const
{
  out.setConstant(std::numeric_limits<double>::infinity());
}

void S2::getTypicalMagnitude_(RefVec out) const { out = typicalMagnitude_; }

void S2::setTypicalMagnitude(double magnitude)
{
  setTypicalMagnitude(Eigen::VectorXd::Constant(tangentDim(), magnitude));
}

void S2::setTypicalMagnitude(const ConstRefVec& out)
{
  testLock();
  typicalMagnitude_ = out;
}

void S2::getTrustMagnitude_(RefVec out) const { out = trustMagnitude_; }

void S2::setTrustMagnitude(const double& magnitude)
{
  setTrustMagnitude(Eigen::VectorXd::Constant(tangentDim(), magnitude));
}

void S2::setTrustMagnitude(const ConstRefVec& out)
{
  testLock();
  trustMagnitude_ = out;
}

long S2::getTypeId() const
{
  long typeId = utils::hash::computeHash("S2");
  return typeId;
}

std::shared_ptr<Manifold> S2::getNewCopy_() const
{
  std::shared_ptr<S2> copy(new S2(*this));
  return copy;
}

bool S2::isSameTopology(const Manifold& other) const
{
  if (dynamic_cast<const S2*>(&other))
    return true;
  else
    return false;
}
}
