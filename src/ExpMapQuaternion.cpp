#include <iostream>
#include <boost/math/special_functions/sinc.hpp>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <manifolds/defs.h>
#include <manifolds/ExpMapQuaternion.h>
#include <manifolds/pgs_assert.h>

namespace pgs
{
  const double ExpMapQuaternion::prec = 1e-8; //TODO Should be sqrt(sqrt(machine precision))

  void ExpMapQuaternion::plus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v)
  {
    DisplayType q;
    exponential(q,v);
    Eigen::Map<DisplayType>(out.data()) = (Eigen::Map<const DisplayType>(x.data()))*q;
  }

  void ExpMapQuaternion::exponential(DisplayType& q, const ConstRefVec& v)
  {
    pgs_assert(v.size() == 3 && "Increment for expMap must be of size 3");
    double n = v.squaredNorm(); // (theta)^2 (in Grassia)
    pgs_assert(sqrt(n) < M_PI && "Increment for expMap must be of norm at most pi");
    double c; // cos(theta/2) in Grassia
    double s; // sin(theta/2)/theta in Grassia
    if (n < prec)
    {
      c = 1 - n / 8;
      s = 1/2 + n / 48;
    }
    else
    {
      double t = sqrt(n); // theta (in Grassia)
      c = cos(t/2);
      s = sin(t/2) / t;
    }
    q.w() = c;
    q.vec() = s*v;
  }

  void ExpMapQuaternion::minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& y)
  {
    typedef Eigen::Map<const Eigen::Matrix3d> ConstMapMat3;
    DisplayType R(((ConstMapMat3(y.data())).transpose())*(ConstMapMat3(x.data())));
    logarithm(out,R);
  }

  void ExpMapQuaternion::invMap_(RefVec out, const ConstRefVec& x)
  {
    typedef Eigen::Map<const Eigen::Matrix3d> ConstMapMat3;
    DisplayType R(ConstMapMat3(x.data()));
    logarithm(out,R);
  }

  void ExpMapQuaternion::logarithm(RefVec out, const DisplayType& R)
  {
    Eigen::Vector3d v(R.vec());
    out = v;
  }

  void ExpMapQuaternion::setIdentity_(RefVec out)
  {
    out << 1,0,0,0;
  }

  bool ExpMapQuaternion::isValidInit_(const Eigen::VectorXd& val)
  {
    typedef Eigen::Map<const Eigen::Quaterniond> toQuat;
    bool out(val.size()==4);
    double norm = toQuat(val.data()).norm();
    out = out && (fabs(norm - 1) < prec);
    return out;
  }

  Eigen::Matrix<double, 4, 3> ExpMapQuaternion::diffMap_(const ConstRefVec& )
  {
    Eigen::Matrix<double, 4, 3> J;
    J.setZero();
    return J;
  }

  void ExpMapQuaternion::applyDiffMap_(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x, ReusableTemporaryMap& m)
  {
    pgs_assert(in.cols() == OutputDim_ && "Dimensions mismatch" );
    Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> a = m.getMap(in.rows(),3);
    a.noalias() = in*diffMap_(x);
    out = a;
  }

  Eigen::Matrix<double, 3, 4> ExpMapQuaternion::diffInvMap_(const ConstRefVec& )
  {
    Eigen::Matrix<double, 3, 4> J;
    J.setZero();
    return J;
  }

  void ExpMapQuaternion::applyDiffInvMap_(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x, ReusableTemporaryMap& m)
  {
    pgs_assert(in.cols() == InputDim_ && "Dimensions mismatch" );
    Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> a = m.getMap(in.rows(),OutputDim_);
    a.noalias() = in*diffInvMap_(x);
    out = a;
  }

  void ExpMapQuaternion::applyTransport_(RefMat , const ConstRefMat&, const ConstRefVec&, const ConstRefVec& , ReusableTemporaryMap& )
  {
  }

  void ExpMapQuaternion::applyInvTransport_(RefMat , const ConstRefMat& , const ConstRefVec&, const ConstRefVec& , ReusableTemporaryMap& )
  {
  }

  void ExpMapQuaternion::tangentConstraint_(RefMat, const ConstRefVec&)
  {
    //out is 0xt, no need to fill it
  }

  bool ExpMapQuaternion::isInTxM_(const ConstRefVec&, const ConstRefVec&)
  {
    return true;
  }

  void ExpMapQuaternion::forceOnTxM_(RefVec out, const ConstRefVec& in, const ConstRefVec&)
  {
    out = in;
  }

}


