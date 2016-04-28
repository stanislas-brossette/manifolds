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
#include <Eigen/Geometry>
#include <manifolds/defs.h>
#include <manifolds/ExpMapQuaternion.h>
#include <manifolds/mnf_assert.h>

namespace mnf
{
namespace utils
{
ReverseQuaternion::ReverseQuaternion(double* data)
    : Eigen::Quaterniond(data[0], data[1], data[2], data[3]), refData_(data)
{
}

void ReverseQuaternion::print() const
{
  std::cout << "w: " << w() << std::endl;
  std::cout << "x: " << x() << std::endl;
  std::cout << "y: " << y() << std::endl;
  std::cout << "z: " << z() << std::endl;
}

void ReverseQuaternion::writeChanges()
{
  refData_[0] = w();
  refData_[1] = x();
  refData_[2] = y();
  refData_[3] = z();
}

ReverseQuaternion::~ReverseQuaternion() { writeChanges(); }

ReverseQuaternion& ReverseQuaternion::operator=(const Eigen::Quaterniond& quat)
{
  w() = quat.w();
  vec() = quat.vec();
  return *this;
}

ConstReverseQuaternion::ConstReverseQuaternion(const double* data)
    : Eigen::Quaterniond(data[0], data[1], data[2], data[3])
{
}
void ConstReverseQuaternion::print() const
{
  std::cout << "w: " << w() << std::endl;
  std::cout << "x: " << x() << std::endl;
  std::cout << "y: " << y() << std::endl;
  std::cout << "z: " << z() << std::endl;
}
}

typedef utils::ReverseQuaternion toQuat;
typedef utils::ConstReverseQuaternion toConstQuat;
// typedef Eigen::Map< Eigen::Quaterniond > toQuat;
// typedef Eigen::Map< const Eigen::Quaterniond > toConstQuat;
const double ExpMapQuaternion::prec =
    1e-8;  // TODO Should be sqrt(sqrt(machine precision))
#if defined(_MSC_FULL_VER) && _MSC_VER <= 1900
char ExpMapQuaternion::hashName[] = "ExpMapQuaternion";
#endif

void ExpMapQuaternion::retractation_(RefVec out, const ConstRefVec& x,
                                     const ConstRefVec& v)
{
  OutputType q;
  exponential(q, v);
  toQuat(out.data()) =
      (toConstQuat(x.data())) * (toConstQuat(q.data()));  // out = x*exp(v)
}

void ExpMapQuaternion::exponential(OutputType& q, const ConstRefVec& v)
{
  mnf_assert(v.size() == 3 && "Increment for expMap must be of size 3");
  double n2 = v.squaredNorm();  // (theta)^2 (in Grassia)
  mnf_assert(sqrt(n2) < M_PI &&
             "Increment for expMap must be of norm at most pi");
  double s;  // sin(theta/2)/theta in Grassia
  if (n2 < prec)
  {
    toQuat(q.data()).w() =
        1 + (-1 + n2 / 48) * (n2 / 8);  // cos(theta/2) in Grassia
    s = (1 + (-1 + 0.0125 * n2) * n2 / 24) / 2;
  }
  else
  {
    double t = sqrt(n2);  // theta (in Grassia)
    toQuat(q.data()).w() = cos(0.5 * t);
    s = sin(0.5 * t) / t;
  }
  toQuat(q.data()).vec() = s * v;
}

void ExpMapQuaternion::pseudoLog_(RefVec out, const ConstRefVec& x,
                                  const ConstRefVec& y)
{
  Eigen::Vector4d tmp;
  toQuat q(tmp.data());
  const toConstQuat xQ(x.data());
  const toConstQuat yQ(y.data());
  q = xQ.inverse() * yQ;  // TODO double-check that formula
  q.writeChanges();
  logarithm(out, tmp);
}

void ExpMapQuaternion::pseudoLog0_(RefVec out, const ConstRefVec& x)
{
  logarithm(out, x);
}

Eigen::Vector3d ExpMapQuaternion::getLog(const ConstRefVec& x, const ConstRefVec& y, Eigen::Vector4d& tmp)
{
  toQuat q(tmp.data());
  const toConstQuat xQ(x.data());
  const toConstQuat yQ(y.data());
  q = xQ.inverse() * yQ;  // TODO double-check that formula
  q.writeChanges();
  Eigen::Vector3d d;
  logarithm(d, tmp);
  return d;
}

double ExpMapQuaternion::distance_(const ConstRefVec& x,
                              const ConstRefVec& y)
{
  Eigen::Vector4d tmp;
  return getLog(x, y, tmp).norm();
}

double ExpMapQuaternion::squaredDistance_(const ConstRefVec& x,
                              const ConstRefVec& y)
{
  Eigen::Vector4d tmp;
  return getLog(x, y, tmp).squaredNorm();
}

double ExpMapQuaternion::squaredDistance_(const ConstRefVec& x,
                              const ConstRefVec& y, const ConstRefVec& w)
{
  Eigen::Vector4d tmp;
  Eigen::Vector3d d = getLog(x,y, tmp);
  d = (w.array() * d.array()).matrix();
  return d.squaredNorm();
}

Eigen::Matrix<double, 4, 4> diffInvXTimesYdx(const ConstRefVec& x,
                                             const ConstRefVec& y)
{
  Eigen::Matrix<double, 4, 4> diffXtYdx;

  toConstQuat xQ(x.data());
  toConstQuat yQ(y.data());
  double nxQ = xQ.w()*xQ.w()+xQ.x()*xQ.x()+xQ.y()*xQ.y()+xQ.z()*xQ.z();
  double nxQ2 = pow(nxQ,2.0);

  diffXtYdx(0, 0) = 1.0/nxQ2*(-xQ.w()*(xQ.x()*yQ.x()*2.0+xQ.y()*yQ.y()*2.0+xQ.z()*yQ.z()*2.0)-(xQ.w()*xQ.w())*yQ.w()+yQ.w()*(xQ.x()*xQ.x())+yQ.w()*(xQ.y()*xQ.y())+yQ.w()*(xQ.z()*xQ.z()));
  diffXtYdx(0, 1) = 1.0/nxQ2*(-xQ.x()*(xQ.w()*yQ.w()*2.0+xQ.y()*yQ.y()*2.0+xQ.z()*yQ.z()*2.0)+(xQ.w()*xQ.w())*yQ.x()-(xQ.x()*xQ.x())*yQ.x()+yQ.x()*(xQ.y()*xQ.y())+yQ.x()*(xQ.z()*xQ.z()));
  diffXtYdx(0, 2) = 1.0/nxQ2*(-xQ.y()*(xQ.w()*yQ.w()*2.0+xQ.x()*yQ.x()*2.0+xQ.z()*yQ.z()*2.0)+(xQ.w()*xQ.w())*yQ.y()+(xQ.x()*xQ.x())*yQ.y()-(xQ.y()*xQ.y())*yQ.y()+yQ.y()*(xQ.z()*xQ.z()));
  diffXtYdx(0, 3) = 1.0/nxQ2*(-xQ.z()*(xQ.w()*yQ.w()*2.0+xQ.x()*yQ.x()*2.0+xQ.y()*yQ.y()*2.0)+(xQ.w()*xQ.w())*yQ.z()+(xQ.x()*xQ.x())*yQ.z()+(xQ.y()*xQ.y())*yQ.z()-(xQ.z()*xQ.z())*yQ.z());
  diffXtYdx(1, 0) = 1.0/nxQ2*(xQ.w()*(yQ.w()*xQ.x()*2.0+xQ.y()*yQ.z()*2.0-yQ.y()*xQ.z()*2.0)-(xQ.w()*xQ.w())*yQ.x()+(xQ.x()*xQ.x())*yQ.x()+yQ.x()*(xQ.y()*xQ.y())+yQ.x()*(xQ.z()*xQ.z()));
  diffXtYdx(1, 1) = -1.0/nxQ2*(xQ.x()*(xQ.w()*yQ.x()*2.0-xQ.y()*yQ.z()*2.0+yQ.y()*xQ.z()*2.0)+(xQ.w()*xQ.w())*yQ.w()-yQ.w()*(xQ.x()*xQ.x())+yQ.w()*(xQ.y()*xQ.y())+yQ.w()*(xQ.z()*xQ.z()));
  diffXtYdx(1, 2) = -1.0/nxQ2*(xQ.y()*(xQ.w()*yQ.x()*2.0-yQ.w()*xQ.x()*2.0+yQ.y()*xQ.z()*2.0)+(xQ.w()*xQ.w())*yQ.z()+(xQ.x()*xQ.x())*yQ.z()-(xQ.y()*xQ.y())*yQ.z()+(xQ.z()*xQ.z())*yQ.z());
  diffXtYdx(1, 3) = 1.0/nxQ2*(xQ.z()*(xQ.w()*yQ.x()*-2.0+yQ.w()*xQ.x()*2.0+xQ.y()*yQ.z()*2.0)+(xQ.w()*xQ.w())*yQ.y()+(xQ.x()*xQ.x())*yQ.y()+(xQ.y()*xQ.y())*yQ.y()-yQ.y()*(xQ.z()*xQ.z()));
  diffXtYdx(2, 0) = 1.0/nxQ2*(xQ.w()*(yQ.w()*xQ.y()*2.0-xQ.x()*yQ.z()*2.0+yQ.x()*xQ.z()*2.0)-(xQ.w()*xQ.w())*yQ.y()+(xQ.x()*xQ.x())*yQ.y()+(xQ.y()*xQ.y())*yQ.y()+yQ.y()*(xQ.z()*xQ.z()));
  diffXtYdx(2, 1) = 1.0/nxQ2*(xQ.x()*(xQ.w()*yQ.y()*-2.0+yQ.w()*xQ.y()*2.0+yQ.x()*xQ.z()*2.0)+(xQ.w()*xQ.w())*yQ.z()-(xQ.x()*xQ.x())*yQ.z()+(xQ.y()*xQ.y())*yQ.z()+(xQ.z()*xQ.z())*yQ.z());
  diffXtYdx(2, 2) = -1.0/nxQ2*(xQ.y()*(xQ.w()*yQ.y()*2.0+xQ.x()*yQ.z()*2.0-yQ.x()*xQ.z()*2.0)+(xQ.w()*xQ.w())*yQ.w()+yQ.w()*(xQ.x()*xQ.x())-yQ.w()*(xQ.y()*xQ.y())+yQ.w()*(xQ.z()*xQ.z()));
  diffXtYdx(2, 3) = -1.0/nxQ2*(xQ.z()*(xQ.w()*yQ.y()*2.0-yQ.w()*xQ.y()*2.0+xQ.x()*yQ.z()*2.0)+(xQ.w()*xQ.w())*yQ.x()+(xQ.x()*xQ.x())*yQ.x()+yQ.x()*(xQ.y()*xQ.y())-yQ.x()*(xQ.z()*xQ.z()));
  diffXtYdx(3, 0) = 1.0/nxQ2*(xQ.w()*(yQ.w()*xQ.z()*2.0+xQ.x()*yQ.y()*2.0-yQ.x()*xQ.y()*2.0)-(xQ.w()*xQ.w())*yQ.z()+(xQ.x()*xQ.x())*yQ.z()+(xQ.y()*xQ.y())*yQ.z()+(xQ.z()*xQ.z())*yQ.z());
  diffXtYdx(3, 1) = -1.0/nxQ2*(xQ.x()*(xQ.w()*yQ.z()*2.0-yQ.w()*xQ.z()*2.0+yQ.x()*xQ.y()*2.0)+(xQ.w()*xQ.w())*yQ.y()-(xQ.x()*xQ.x())*yQ.y()+(xQ.y()*xQ.y())*yQ.y()+yQ.y()*(xQ.z()*xQ.z()));
  diffXtYdx(3, 2) = 1.0/nxQ2*(xQ.y()*(xQ.w()*yQ.z()*-2.0+yQ.w()*xQ.z()*2.0+xQ.x()*yQ.y()*2.0)+(xQ.w()*xQ.w())*yQ.x()+(xQ.x()*xQ.x())*yQ.x()-yQ.x()*(xQ.y()*xQ.y())+yQ.x()*(xQ.z()*xQ.z()));
  diffXtYdx(3, 3) = -1.0/nxQ2*(xQ.z()*(xQ.w()*yQ.z()*2.0-xQ.x()*yQ.y()*2.0+yQ.x()*xQ.y()*2.0)+(xQ.w()*xQ.w())*yQ.w()+yQ.w()*(xQ.x()*xQ.x())+yQ.w()*(xQ.y()*xQ.y())-yQ.w()*(xQ.z()*xQ.z()));

  return diffXtYdx;
}

Eigen::Matrix<double, 4, 4> diffInvXTimesYdy(const ConstRefVec& x,
                                             const ConstRefVec&)
{
  Eigen::Matrix<double, 4, 4> diffXtYdy;

  toConstQuat xQ(x.data());

  double nxQ = xQ.w()*xQ.w()+xQ.x()*xQ.x()+xQ.y()*xQ.y()+xQ.z()*xQ.z();

  diffXtYdy(0, 0) = xQ.w()/nxQ;
  diffXtYdy(0, 1) = xQ.x()/nxQ;
  diffXtYdy(0, 2) = xQ.y()/nxQ;
  diffXtYdy(0, 3) = xQ.z()/nxQ;
  diffXtYdy(1, 0) = -xQ.x()/nxQ;
  diffXtYdy(1, 1) = xQ.w()/nxQ;
  diffXtYdy(1, 2) = xQ.z()/nxQ;
  diffXtYdy(1, 3) = -xQ.y()/nxQ;
  diffXtYdy(2, 0) = -xQ.y()/nxQ;
  diffXtYdy(2, 1) = -xQ.z()/nxQ;
  diffXtYdy(2, 2) = xQ.w()/nxQ;
  diffXtYdy(2, 3) = xQ.x()/nxQ;
  diffXtYdy(3, 0) = -xQ.z()/nxQ;
  diffXtYdy(3, 1) = xQ.y()/nxQ;
  diffXtYdy(3, 2) = -xQ.x()/nxQ;
  diffXtYdy(3, 3) = xQ.w()/nxQ;

  return diffXtYdy;
}

Eigen::Matrix<double, 1, 4> ExpMapQuaternion::derivDistanceX_(
    const ConstRefVec& x, const ConstRefVec& y)
{
  Eigen::Vector4d tmp;
  Eigen::Vector3d d = getLog(x,y, tmp);
  Eigen::Matrix<double, 1, 3> n = d.transpose() / d.norm();

  Eigen::Matrix<double, 3, 4> diffLog = diffPseudoLog0_(tmp);

  Eigen::Matrix<double, 1, 4> J;
  J = n * diffLog * diffInvXTimesYdx(x, y);
  return J;
}
Eigen::Matrix<double, 1, 4> ExpMapQuaternion::derivDistanceY_(
    const ConstRefVec& x, const ConstRefVec& y)
{
  Eigen::Vector4d tmp;
  Eigen::Vector3d d = getLog(x,y, tmp);
  Eigen::Matrix<double, 1, 3> n = d.transpose() / d.norm();

  Eigen::Matrix<double, 3, 4> diffLog = diffPseudoLog0_(tmp);

  Eigen::Matrix<double, 1, 4> J;
  J = n * diffLog * diffInvXTimesYdy(x, y);
  return J;
}
Eigen::Matrix<double, 1, 4> ExpMapQuaternion::derivSquaredDistanceX_(
    const ConstRefVec& x, const ConstRefVec& y)
{
  Eigen::Vector4d tmp;
  Eigen::Vector3d d = getLog(x,y, tmp);
  Eigen::Matrix<double, 1, 3> n = 2 * d.transpose();

  Eigen::Matrix<double, 3, 4> diffLog = diffPseudoLog0_(tmp);

  Eigen::Matrix<double, 1, 4> J;
  J = n * diffLog * diffInvXTimesYdx(x, y);
  return J;
}
Eigen::Matrix<double, 1, 4> ExpMapQuaternion::derivSquaredDistanceY_(
    const ConstRefVec& x, const ConstRefVec& y)
{
  Eigen::Vector4d tmp;
  Eigen::Vector3d d = getLog(x,y, tmp);
  Eigen::Matrix<double, 1, 3> n = 2 * d.transpose();

  Eigen::Matrix<double, 3, 4> diffLog = diffPseudoLog0_(tmp);

  Eigen::Matrix<double, 1, 4> J;
  J = n * diffLog * diffInvXTimesYdy(x, y);
  return J;
}
Eigen::Matrix<double, 1, 4> ExpMapQuaternion::derivSquaredDistanceX_(
    const ConstRefVec& x, const ConstRefVec& y, const ConstRefVec& w)
{
  Eigen::Vector4d tmp;
  Eigen::Vector3d d = getLog(x,y, tmp);
  d = (w.array() * d.array()).matrix();
  Eigen::Matrix<double, 1, 3> n = 2 * d.transpose();

  Eigen::Matrix<double, 3, 4> diffLog = diffPseudoLog0_(tmp);

  Eigen::Matrix<double, 1, 4> J;
  J = n * diffLog * diffInvXTimesYdx(x, y);
  return J;
}
Eigen::Matrix<double, 1, 4> ExpMapQuaternion::derivSquaredDistanceY_(
    const ConstRefVec& x, const ConstRefVec& y, const ConstRefVec& w)
{
  Eigen::Vector4d tmp;
  Eigen::Vector3d d = getLog(x,y, tmp);
  d = (w.array() * d.array()).matrix();
  Eigen::Matrix<double, 1, 3> n = 2 * d.transpose();

  Eigen::Matrix<double, 3, 4> diffLog = diffPseudoLog0_(tmp);

  Eigen::Matrix<double, 1, 4> J;
  J = n * diffLog * diffInvXTimesYdy(x, y);
  return J;
}

void ExpMapQuaternion::logarithm(RefVec out, const OutputType& v)
{
  const toConstQuat vQ(v.data());
  double n2 = vQ.vec().squaredNorm();
  double n = sqrt(n2);

  if (n < prec && vQ.w() != 0)
    out = (2 / vQ.w()) * vQ.vec();
  else
    out = atan2(2 * n * vQ.w(), vQ.w() * vQ.w() - n2) / n * vQ.vec();
}

void ExpMapQuaternion::setZero_(RefVec out)
{
  toQuat(out.data()).setIdentity();
}

bool ExpMapQuaternion::isInM_(const Eigen::VectorXd& val, double prec)
{
  bool out(val.size() == 4);
  double norm = toConstQuat(val.data()).norm();
  out = out && (fabs(norm - 1) < prec);
  return out;
}

void ExpMapQuaternion::forceOnM_(RefVec out, const ConstRefVec& in)
{
  toConstQuat inQuat(in.data());
  toQuat outQuat(out.data());
  outQuat = inQuat;
  outQuat.normalize();
  outQuat.writeChanges();
}

void ExpMapQuaternion::getIdentityOnTxM_(RefMat out, const ConstRefVec&)
{
  out.setIdentity();
}

Eigen::Matrix<double, 4, 3> ExpMapQuaternion::diffRetractation_(
    const ConstRefVec& x)
{
  toConstQuat xQ(x.data());
  Eigen::Matrix<double, 4, 3> J;
  // This matrix is written in the (w, x, y, z) convention
  // for quaternion notation.
  J << -0.5*xQ.x(), -0.5*xQ.y(), -0.5*xQ.z(),
        0.5*xQ.w(), -0.5*xQ.z(),  0.5*xQ.y(),
        0.5*xQ.z(),  0.5*xQ.w(), -0.5*xQ.x(),
       -0.5*xQ.y(),  0.5*xQ.x(),  0.5*xQ.w();
  return J;
}

void ExpMapQuaternion::applyDiffRetractation_(RefMat out, const ConstRefMat& in,
                                              const ConstRefVec& x,
                                              ReusableTemporaryMap& m)
{
  mnf_assert(in.cols() == OutputDim_ && "Dimensions mismatch");
  Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> a = m.getMap(in.rows(), 3);
  a.noalias() = in * diffRetractation_(x);
  out = a;
}

Eigen::Matrix<double, 3, 4> ExpMapQuaternion::diffPseudoLog0_(
    const ConstRefVec& v)
{
  const toConstQuat vQ(v.data());
  double n2 = vQ.vec().squaredNorm();
  double n = sqrt(n2);
  Eigen::Matrix<double, 3, 4> J;

  if (n < prec && vQ.w() != 0)
  {
    double a = 2 / vQ.w();
    double b = -2 / (vQ.w() * vQ.w());
    // This matrix is written in the (w, x, y, z) convention
    // for quaternion notation.
    J << b*vQ.x(), a, 0, 0,
         b*vQ.y(), 0, a, 0,
         b*vQ.z(), 0, 0, a;
  }
  else
  {
    // log(x,y,z,w) = f(x,y,z,w)*[x;y;z]
    double f = atan2(2 * n * vQ.w(), vQ.w() * vQ.w() - n2) / n;
    // df/dx = (x*atan((2*w*(x^2 + y^2 + z^2)^(1/2))/(- w^2 + x^2 + y^2 +
    // z^2)))/(x^2 + y^2 + z^2)^(3/2) - ((2*w*x)/((x^2 + y^2 + z^2)^(1/2)*(- w^2
    // + x^2 + y^2 + z^2)) - (4*w*x*(x^2 + y^2 + z^2)^(1/2))/(- w^2 + x^2 + y^2
    // + z^2)^2)/(((4*w^2*(x^2 + y^2 + z^2))/(- w^2 + x^2 + y^2 + z^2)^2 +
    // 1)*(x^2 + y^2 + z^2)^(1/2))
    // df/dy = (y*atan((2*w*(x^2 + y^2 + z^2)^(1/2))/(- w^2 + x^2 + y^2 +
    // z^2)))/(x^2 + y^2 + z^2)^(3/2) - ((2*w*y)/((x^2 + y^2 + z^2)^(1/2)*(- w^2
    // + x^2 + y^2 + z^2)) - (4*w*y*(x^2 + y^2 + z^2)^(1/2))/(- w^2 + x^2 + y^2
    // + z^2)^2)/(((4*w^2*(x^2 + y^2 + z^2))/(- w^2 + x^2 + y^2 + z^2)^2 +
    // 1)*(x^2 + y^2 + z^2)^(1/2))
    // df/dz = (z*atan((2*w*(x^2 + y^2 + z^2)^(1/2))/(- w^2 + x^2 + y^2 +
    // z^2)))/(x^2 + y^2 + z^2)^(3/2) - ((2*w*z)/((x^2 + y^2 + z^2)^(1/2)*(- w^2
    // + x^2 + y^2 + z^2)) - (4*w*z*(x^2 + y^2 + z^2)^(1/2))/(- w^2 + x^2 + y^2
    // + z^2)^2)/(((4*w^2*(x^2 + y^2 + z^2))/(- w^2 + x^2 + y^2 + z^2)^2 +
    // 1)*(x^2 + y^2 + z^2)^(1/2))

    // g = (atan((2*w*n)/(- w^2 + n2)))/(n2)^(3/2) - ((2*w)/(n*(- w^2 + n2)) -
    // (4*w*n)/(- w^2 + n2)^2)/(((4*w^2*n2)/(- w^2 + n2)^2 + 1)*n)
    // g = (atan((2*w*n)/(- w^2 + n2)))/(n2*n) - ((2*w/n2)*(- w^2 + n2) -
    // 4*w)/(4*w^2*n2+(- w^2 + n2)^2)
    // df/dx = g*x = (x*atan((2*w*n)/(- w^2 + n2)))/(n2)^(3/2) - ((2*w*x)/(n*(-
    // w^2 + n2)) - (4*w*x*n)/(- w^2 + n2)^2)/(((4*w^2*n2)/(- w^2 + n2)^2 +
    // 1)*n)
    // df/dy = g*y = (y*atan((2*w*n)/(- w^2 + n2)))/(n2)^(3/2) - ((2*w*y)/(n*(-
    // w^2 + n2)) - (4*w*y*n)/(- w^2 + n2)^2)/(((4*w^2*n2)/(- w^2 + n2)^2 +
    // 1)*n)
    // df/dz = g*z = (z*atan((2*w*n)/(- w^2 + n2)))/(n2)^(3/2) - ((2*w*z)/(n*(-
    // w^2 + n2)) - (4*w*z*n)/(- w^2 + n2)^2)/(((4*w^2*n2)/(- w^2 + n2)^2 +
    // 1)*n)
    // df/dw = -2/(w²+x²+y²+z²)
    double g = (2 * vQ.w() - f) / n2;
    double dfdw = -2;
    /*
     * J = [ g.x²+f, g.y.x, g.z.x, df/dw.x]
     *     [ g.x.y, g.y²+f, g.z.y, df/dw.y]
     *     [ g.x.z, g.y.z, g.z²+f, df/dw.z]
    */
    // This matrix is written in the (w, x, y, z) convention
    // for quaternion notation.
    // J << dfdw*vQ.x(), g*vQ.x()*vQ.x()+f, g*vQ.y()*vQ.x(), g*vQ.z()*vQ.x(),
    // 	   dfdw*vQ.y(), g*vQ.x()*vQ.y(), g*vQ.y()*vQ.y()+f, g*vQ.z()*vQ.y(),
    // 	   dfdw*vQ.z(), g*vQ.x()*vQ.z(), g*vQ.y()*vQ.z(), g*vQ.z()*vQ.z()+f;

    J.col(0) = dfdw * vQ.vec();
    J.rightCols<3>() = f * Eigen::Matrix3d::Identity();
    J.rightCols<3>().noalias() += g * vQ.vec() * vQ.vec().transpose();
  }
  return J;
}

void ExpMapQuaternion::applyDiffPseudoLog0_(RefMat out, const ConstRefMat& in,
                                            const ConstRefVec& x,
                                            ReusableTemporaryMap& m)
{
  mnf_assert(in.cols() == InputDim_ && "Dimensions mismatch");
  Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> a =
      m.getMap(in.rows(), OutputDim_);
  a.noalias() = in * diffPseudoLog0_(x);
  out = a;
}

// void ExpMapQuaternion::applyTransport_(RefMat out, const ConstRefMat& in,
// const ConstRefVec&, const ConstRefVec& v, ReusableTemporaryMap& m)
void ExpMapQuaternion::applyTransport_(RefMat out, const ConstRefMat& in,
                                       const ConstRefVec&, const ConstRefVec&)
{
  // OutputType E;
  // exponential(E,v);
  // Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> a = m.getMap(InputDim_,
  // in.cols());
  // a.noalias() = E*in;
  // out = a;
  // TODO Make sure that out=in is OK here...
  out = in;
}

// void ExpMapQuaternion::applyInvTransport_(RefMat out, const ConstRefMat& in,
// const ConstRefVec&, const ConstRefVec& v, ReusableTemporaryMap& m)
void ExpMapQuaternion::applyInvTransport_(RefMat out, const ConstRefMat& in,
                                          const ConstRefVec&,
                                          const ConstRefVec&)
{
  // OutputType E;
  // exponential(E,v);
  // Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> a = m.getMap(in.rows(),
  // InputDim_);
  // a.noalias() = in*(E.transpose());
  // out = a;
  // TODO Make sure that out=in is OK here...
  out = in;
}

void ExpMapQuaternion::applyInvTransportOnTheRight_(RefMat out,
                                                    const ConstRefMat& in,
                                                    const ConstRefVec&,
                                                    const ConstRefVec&)
{
  out = in;
}

void ExpMapQuaternion::tangentConstraint_(RefMat, const ConstRefVec&)
{
  // out is 0xt, no need to fill it
}

bool ExpMapQuaternion::isInTxM_(const ConstRefVec&, const ConstRefVec&,
                                const double&)
{
  return true;
}

void ExpMapQuaternion::forceOnTxM_(RefVec out, const ConstRefVec& in,
                                   const ConstRefVec&)
{
  out = in;
}
}
