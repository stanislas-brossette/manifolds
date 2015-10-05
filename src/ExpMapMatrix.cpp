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

#include <iostream>
#include <boost/math/special_functions/sinc.hpp>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <manifolds/defs.h>
#include <manifolds/ExpMapMatrix.h>
#include <manifolds/mnf_assert.h>

namespace utility
{
// areOverlappingData tests if the data pointed by a and b are overlapping, in
// which case, some aliasing could appear.
// TODO, explain more.
// This function assumes that there wasn't any copies of a or b before the
// call.
bool areOverlappingData(const mnf::ConstRefVec& a, const mnf::ConstRefVec& b)
{
  bool res = false;
  if (&a.coeff(a.rows() - 1) < &b.coeff(0) ||
      &b.coeff(b.rows() - 1) < &a.coeff(0))
  {
    return false;
  }

  for (int i = 0; i < a.rows(); ++i)
  {
    for (int j = i; j < b.rows(); ++j)
    {
      res = res || (&a.coeff(i) == &b.coeff(j));
    }
  }
  return res;
}
}

namespace mnf
{
const double ExpMapMatrix::prec = 1e-8;
typedef Eigen::Map<const Eigen::Matrix3d> toConstMat3;
typedef Eigen::Map<Eigen::Matrix3d> toMat3;
#if defined(_MSC_FULL_VER) && _MSC_VER < 1900
char ExpMapMatrix::hashName[] = "ExpMapMatrix";
#endif

void ExpMapMatrix::retractation_(RefVec out, const ConstRefVec& x,
                                 const ConstRefVec& v)
{
  OutputType E;
  exponential(E, v);
  toMat3(out.data()) = (toConstMat3(x.data())) * E;
}

void ExpMapMatrix::exponential(OutputType& E, const ConstRefVec& v)
{
  mnf_assert(v.size() == 3 && "Increment for expMap must be of size 3");
  double n = v.squaredNorm();
  mnf_assert(sqrt(n) < M_PI &&
             "Increment for expMap must be of norm at most pi");
  double c, s;
  if (n < prec)
  {
    c = 0.5 - n / 24;
    s = 1 - n / 6;  // TODO Check this value. Because there is a mistake here in
                    // Grassia
  }
  else
  {
    double t = sqrt(n);
    c = (1 - cos(t)) / n;
    s = sin(t) / t;
  }
  E << 1 - c*(v.y() * v.y() + v.z() * v.z()), -s * v.z() + c * v.x() * v.y(),
      s * v.y() + c * v.x() * v.z(), s * v.z() + c * v.x() * v.y(),
      1 - c * (v.x() * v.x() + v.z() * v.z()), -s * v.x() + c * v.y() * v.z(),
      -s * v.y() + c * v.x() * v.z(), s * v.x() + c * v.y() * v.z(),
      1 - c * (v.x() * v.x() + v.y() * v.y());
}

void ExpMapMatrix::pseudoLog_(RefVec out, const ConstRefVec& x,
                              const ConstRefVec& y)
{
  OutputType R(((toConstMat3(x.data())).transpose()) * (toConstMat3(y.data())));
  logarithm(out, R);
}

void ExpMapMatrix::pseudoLog0_(RefVec out, const ConstRefVec& x)
{
  OutputType R(toConstMat3(x.data()));
  logarithm(out, R);
}

void ExpMapMatrix::logarithm(RefVec out, const OutputType& R)
{
  Eigen::Vector3d v(-R(1, 2), R(0, 2), -R(0, 1));
  double acosTr = std::acos((R.trace() - 1) / 2);
  if (v.norm() < prec)
    out = v;
  else
  {
    OutputType diff(R - R.transpose());
    double coeff = acosTr / (2 * std::sin(acosTr));
    v(0) = diff(2, 1) * coeff;
    v(1) = diff(0, 2) * coeff;
    v(2) = diff(1, 0) * coeff;
    out = v;
  }
}

void ExpMapMatrix::setZero_(RefVec out)
{
  toMat3(out.data()) = Eigen::Matrix3d::Identity();
}

bool ExpMapMatrix::isInM_(const Eigen::VectorXd& val, double)
{
  bool out(val.size() == 9);
  toConstMat3 valMat(val.data());
  double det = valMat.determinant();
  out = out && (fabs(det - 1) < prec);
  out = out && ((valMat.transpose()) * valMat).isIdentity(prec);
  return out;
}

void ExpMapMatrix::forceOnM_(RefVec out, const ConstRefVec& in)
{
  toConstMat3 inMat(in.data());
  mnf_assert(
      (inMat.transpose() * inMat).isApprox(Eigen::Matrix3d::Identity(), 0.1) &&
      fabs(inMat.determinant() - 1.0) < 0.1 &&
      "Provided matrix is too far from being a rotation matrix. You should use "
      "createRandomPoint instead of forceOnM");
  Eigen::Matrix3d A(
      inMat);  // TODO it is bad to copy it here. Maybe in shouldn't be const???
  toMat3 Q(out.data());
  Eigen::Matrix3d R;
  // the following scope contains the code to compute a QR on R3x3
  {
    // k = 0
    R(0, 0) = A.col(0).norm();
    Q.col(0) = A.col(0) / R(0, 0);
    // j = 1
    R(0, 1) = Q.col(0).dot(A.col(1));
    A.col(1) = A.col(1) - Q.col(0) * R(0, 1);
    // j = 2
    R(0, 2) = Q.col(0).dot(A.col(2));
    A.col(2) = A.col(2) - Q.col(0) * R(0, 2);
    // k = 1
    R(1, 1) = A.col(1).norm();
    Q.col(1) = A.col(1) / R(1, 1);
    // j = 2
    R(1, 2) = Q.col(1).dot(A.col(2));
    A.col(2) = A.col(2) - Q.col(1) * R(1, 2);
    // k = 2
    R(2, 2) = A.col(2).norm();
    Q.col(2) = A.col(2) / R(2, 2);
  }
}

void ExpMapMatrix::getIdentityOnTxM_(RefMat out, const ConstRefVec&)
{
  out.setIdentity();
}

Eigen::Matrix<double, 9, 3> ExpMapMatrix::diffRetractation_(
    const ConstRefVec& x)
{
  Eigen::Matrix<double, 9, 3> J;
  J << 0, -x(6), x(3), 0, -x(7), x(4), 0, -x(8), x(5), x(6), 0, -x(0), x(7), 0,
      -x(1), x(8), 0, -x(2), -x(3), x(0), 0, -x(4), x(1), 0, -x(5), x(2), 0;
  return J;
}

void ExpMapMatrix::applyDiffRetractation_(RefMat out, const ConstRefMat& in,
                                          const ConstRefVec& x,
                                          ReusableTemporaryMap& m)
{
  mnf_assert(in.cols() == OutputDim_ && "Dimensions mismatch");
  Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> a = m.getMap(in.rows(), 3);
  a.noalias() = in * diffRetractation_(x);
  out = a;
}

Eigen::Matrix<double, 3, 9> ExpMapMatrix::diffPseudoLog0_(const ConstRefVec& R)
{
  Eigen::Matrix<double, 3, 9> J;
  J.setZero();
  // Valid approximation of the log when v<<1
  Eigen::Vector3d v((R(5) - R(7)) / 2, (R(6) - R(2)) / 2, (R(1) - R(3)) / 2);

  double trR = R(0) + R(4) + R(8);  // Trace of R;
  double trRm1o2 = (trR - 1) / 2;
  double s2 = 1 - trRm1o2 * trRm1o2;
  double f = 1 / boost::math::sinc_pi(
                     acos(trRm1o2));  //=acos(trRm1o2)/sin(acos(trRm1o2));
  double hf = f / 2;
  double df;
  if (v.squaredNorm() < 1e-10)  // Probably should not use this switch value
                                // here
  {
    // TODO: I do not know how to test that. Verification with finite
    // differences proved unefficient. Need to get back to it later.
    //
    // Here we use the Taylor approximation of the diff of the log
    double x = trR - 3;
    df = trR / 15 - (3 * (x * x)) / 140 + (2 * (x * x * x)) / 315 -
         (5 * (x * x * x * x)) / 2772 - 11 / 30;
  }
  else
  {
    df = (trRm1o2 * f - 1) / (2 * s2);
  }
  J.col(0) = df * v;
  J.col(1) = Eigen::Vector3d(0, 0, hf);
  J.col(2) = Eigen::Vector3d(0, -hf, 0);
  J.col(3) = Eigen::Vector3d(0, 0, -hf);
  J.col(4) = J.col(0);
  J.col(5) = Eigen::Vector3d(hf, 0, 0);
  J.col(6) = Eigen::Vector3d(0, hf, 0);
  J.col(7) = Eigen::Vector3d(-hf, 0, 0);
  J.col(8) = J.col(0);
  return J;
}

void ExpMapMatrix::applyDiffPseudoLog0_(RefMat out, const ConstRefMat& in,
                                        const ConstRefVec& x,
                                        ReusableTemporaryMap& m)
{
  mnf_assert(in.cols() == InputDim_ && "Dimensions mismatch");
  Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> a =
      m.getMap(in.rows(), OutputDim_);
  a.noalias() = in * diffPseudoLog0_(x);
  out = a;
}

void ExpMapMatrix::applyTransport_(RefMat out, const ConstRefMat& in,
                                   const ConstRefVec&, const ConstRefVec&)
{
  // TODO Make sure that out=in is OK here...
  out = in;
}

void ExpMapMatrix::applyInvTransport_(RefMat out, const ConstRefMat& in,
                                      const ConstRefVec&, const ConstRefVec&)
{
  // TODO Make sure that out=in is OK here...
  out = in;
}

void ExpMapMatrix::applyInvTransportOnTheRight_(RefMat out,
                                                const ConstRefMat& in,
                                                const ConstRefVec&,
                                                const ConstRefVec&)
{
  out = in;
}

void ExpMapMatrix::tangentConstraint_(RefMat, const ConstRefVec&)
{
  // out is 0xt, no need to fill it
}

bool ExpMapMatrix::isInTxM_(const ConstRefVec&, const ConstRefVec&,
                            const double&)
{
  return true;
}

void ExpMapMatrix::forceOnTxM_(RefVec out, const ConstRefVec& in,
                               const ConstRefVec&)
{
  out = in;
}
}

