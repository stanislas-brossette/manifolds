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
#include <manifolds/S2_Base.h>
#include <manifolds/mnf_assert.h>

namespace mnf
{
  S2_Base::S2_Base()
    : Manifold_Base(2, 3, 3)
  {
    name() = "S2";
    setTypicalMagnitude(Eigen::Vector3d::Constant(M_PI));
    setTrustMagnitude(Eigen::Vector3d::Constant(M_PI));
  }

  S2_Base::S2_Base(double magnitude)
    : Manifold_Base(2, 3, 3)
  {
    name() = "S2";
    setTypicalMagnitude(Eigen::Vector3d::Constant( magnitude));
    setTrustMagnitude(Eigen::Vector3d::Constant( magnitude));
  }
  S2_Base::S2_Base(const ConstRefVec& magnitude)
    : Manifold_Base(2, 3, 3)
  {
    mnf_assert(magnitude.size() == 3 && "magnitude on R^n must be of size n");
    name() = "S2";
    setTypicalMagnitude (magnitude);
    setTrustMagnitude (magnitude);
  }

  bool S2_Base::isInM_(const Eigen::VectorXd& val , double prec) const
  {
    bool out(fabs(val.lpNorm<2>() - 1.0) < prec);
    return out;
  }

  void S2_Base::forceOnM_(RefVec out, const ConstRefVec& in) const
  {
    out = in/in.lpNorm<2>();
  }

  void S2_Base::getIdentityOnTxM_(RefMat out, const ConstRefVec& x) const
  {
    out.setIdentity();
    Eigen::Matrix3d xxt;
    xxt.noalias()= x*x.transpose();
    out -= xxt;
  }

  size_t S2_Base::numberOfSubmanifolds() const
  {
    return 1;
  }

  bool S2_Base::isElementary() const
  {
    return true;
  }

  std::shared_ptr<const Manifold_Base> S2_Base::operator()(size_t) const
  {
    return shared_from_this();
  }

  std::string S2_Base::toString(const ConstRefVec& val, const std::string& prefix, const Eigen::IOFormat& fmt) const
  {
    std::stringstream ss;
    ss << prefix << val.transpose().format(fmt);
    return ss.str();
  }

  void S2_Base::createRandomPoint_(RefVec out, double ) const
  {
    rand(out);
  }

  void S2_Base::retractation_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    Eigen::Vector3d sum;
    sum = x+v;
    out = sum/sum.lpNorm<2>();
  }

  void S2_Base::pseudoLog_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    logarithm(out, x, y);
  }

  void S2_Base::logarithm (RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    Eigen::Vector3d diff, projDiff;
    diff = y-x;
    projVec(projDiff, diff, x);
    out = distance(x,y)*projDiff/projDiff.lpNorm<2>();
  }

  double S2_Base::distance (const ConstRefVec& x, const ConstRefVec& y) const
  {
    return acos(x.dot(y));
  }

  void S2_Base::projRows(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    for (Index i = 0; i < in.rows(); ++i)
      out.row(i) = in.row(i) - (x*in.row(i)).trace()*x.transpose();
  }
  void S2_Base::projCols(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    for (Index i = 0; i < in.cols(); ++i)
      projVec(out.col(i), in.col(i), x);
  }
  void S2_Base::projVec(RefVec out, const ConstRefVec& in, const ConstRefVec& x) const
  {
      out = in - (x.dot(in))*x;
  }

  void S2_Base::rand(RefVec out) const
  {
    mnf_assert( out.size() == 3 && "Wrong size of output for S2::rand");
    out = Eigen::Vector3d::Random();
    double nout = out.lpNorm<2>();
    out = out/nout;
  }

  void S2_Base::randVec(RefVec out, const ConstRefVec& x) const
  {
    mnf_assert( out.size() == 3 && "Wrong size of output for S2::rand");
    mnf_assert( x.size() == 3 && "Wrong size of output for S2::rand");
    out = Eigen::Vector3d::Random();
    projVec(out, out, x);
    double nout = out.lpNorm<2>();
    out = out/nout;
  }
  Eigen::Vector3d S2_Base::randVec(const ConstRefVec& x) const
  {
    Eigen::Vector3d out;
    randVec(out, x);
    return out;
  }

  void S2_Base::pseudoLog0_(RefVec , const ConstRefVec& ) const
  {
    throw std::runtime_error("Unimplemented method in S2");
  }

  void S2_Base::setZero_(RefVec out) const
  {
    std::cout << "Warning, you are using SetZero on S2. Returning an arbitrary point on S2, which is (1, 0, 0)" << std::endl;
    out << 1.0, 0.0, 0.0;
  }

  Eigen::MatrixXd S2_Base::diffRetractation_(const ConstRefVec& out) const
  {
    throw std::runtime_error("Unimplemented method in S2");
    return out;
  }

  void S2_Base::applyDiffRetractation_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    projRows(out, in, x);
  }

  Eigen::MatrixXd S2_Base::diffPseudoLog0_(const ConstRefVec& out) const
  {
    throw std::runtime_error("Unimplemented method in S2");
    return out;
  }

  void S2_Base::applyDiffPseudoLog0_(RefMat , const ConstRefMat& , const ConstRefVec& ) const
  {
    throw std::runtime_error("Unimplemented method in S2");
  }

  void S2_Base::applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    Eigen::Vector3d y;
    retractation(y, x, v);
    Eigen::Matrix3d R = utils::computeRotBetweenVec(x,y);
    Eigen::Vector3d tmp;
    for (Index i = 0; i < in.cols(); ++i)
    {
      tmp = R*in.col(i);
      out.col(i) = tmp;
    }
  }

  void S2_Base::applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    Eigen::Vector3d y;
    retractation(y, x, v);
    Eigen::Matrix3d R = utils::computeRotBetweenVec(x,y);
    Eigen::Vector3d tmp;
    for (Index i = 0; i < in.cols(); ++i)
    {
      tmp = R.transpose()*in.col(i);
      out.col(i) = tmp;
    }
  }

  void S2_Base::applyInvTransportOnTheRight_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    Eigen::Vector3d y;
    retractation(y, x, v);
    Eigen::Matrix3d R = utils::computeRotBetweenVec(x,y);
    Eigen::Matrix<double, 3, 1> tmp;
    for (Index i = 0; i < in.rows(); ++i)
    {
      tmp = in.row(i)*R.transpose();
      out.row(i) = tmp;
    }
  }

  void S2_Base::tangentConstraint_(RefMat out, const ConstRefVec& x) const
  {
    out = x.transpose();
  }

  bool S2_Base::isInTxM_(const ConstRefVec& x, const ConstRefVec& v, const double& prec) const
  {
    bool out = (fabs(x.dot(v)) < prec);
    return out;
  }

  void S2_Base::forceOnTxM_(RefVec out, const ConstRefVec& v, const ConstRefVec& x) const
  {
    out = v - x.dot(v)*x;
  }

  void S2_Base::limitMap_(RefVec out) const
  {
    out.setConstant(std::numeric_limits<double>::infinity());
  }

  long S2_Base::getTypeId() const
  {
    long typeId = utils::hash::computeHash("S2");
    return typeId;
  }

  Manifold_Base_ptr S2_Base::getNewCopy() const
  {
    S2_Base_ptr copy(new S2_Base(*this));
    copy->instanceId_ = this->instanceId_;

    return copy;
  }
}
