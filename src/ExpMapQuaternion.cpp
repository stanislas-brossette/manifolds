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
#include <manifolds/pgs_assert.h>

namespace pgs
{
  typedef Eigen::Map< Eigen::Quaterniond > toQuat;
  typedef Eigen::Map< const Eigen::Quaterniond > toConstQuat;
  const double ExpMapQuaternion::prec = 1e-8; //TODO Should be sqrt(sqrt(machine precision))

  void ExpMapQuaternion::retractation_(RefVec out, const ConstRefVec& x, const ConstRefVec& v)
  {
    OutputType q;
    exponential(q,v);
    toQuat(out.data()) = (toConstQuat(x.data()))*(toConstQuat(q.data())); //out = x*exp(v)
  }

  void ExpMapQuaternion::exponential(OutputType& q, const ConstRefVec& v)
  {
    pgs_assert(v.size() == 3 && "Increment for expMap must be of size 3");
    double n2 = v.squaredNorm(); // (theta)^2 (in Grassia)
    pgs_assert(sqrt(n2) < M_PI && "Increment for expMap must be of norm at most pi");
    double s; // sin(theta/2)/theta in Grassia
    if (n2 < prec)
    {
      toQuat(q.data()).w() = 1 + (-1 + n2 / 48)*(n2/8);// cos(theta/2) in Grassia
      s = (1+(-1+0.0125*n2)*n2/24)/2;
    }
    else
    {
      double t = sqrt(n2); // theta (in Grassia)
      toQuat(q.data()).w() = cos(0.5*t);
      s = sin(0.5*t) / t;
    }
    toQuat(q.data()).vec() = s*v;
  }

  void ExpMapQuaternion::pseudoLog_(RefVec out, const ConstRefVec& x, const ConstRefVec& y)
  {
    Eigen::Vector4d tmp;
    toQuat q(tmp.data());
    const toConstQuat xQ(x.data());
    const toConstQuat yQ(y.data());
    q = xQ.inverse()*yQ; //TODO double-check that formula
    logarithm(out,tmp);
  }

  void ExpMapQuaternion::pseudoLog0_(RefVec out, const ConstRefVec& x)
  {
    logarithm(out,x);
  }

  void ExpMapQuaternion::logarithm(RefVec out, const OutputType& v)
  {
    const toConstQuat vQ(v.data());
    double n2 = vQ.vec().squaredNorm();
    double n = sqrt(n2);

    if (n < prec && vQ.w()!=0 )
      out = (2/vQ.w())*vQ.vec();
    else
      out = atan2(2 * n * vQ.w(), vQ.w() * vQ.w() - n2) / n * vQ.vec(); 
  }

  void ExpMapQuaternion::setZero_(RefVec out)
  {
    toQuat(out.data()).setIdentity();
  }

  bool ExpMapQuaternion::isInM_(const Eigen::VectorXd& val, const double& )
  {
    bool out(val.size()==4);
    double norm = toConstQuat(val.data()).norm();
    out = out && (fabs(norm - 1) < prec);
    return out;
  }

  void ExpMapQuaternion::forceOnM_(RefVec out, const ConstRefVec& in)
  {
    toConstQuat inQuat(in.data());
    toQuat outQuat(out.data());
    outQuat = inQuat;
    out.normalize();
  }

  Eigen::Matrix<double, 4, 3> ExpMapQuaternion::diffRetractation_(const ConstRefVec& x)
  {
    const Eigen::Map<const Eigen::Quaterniond> xQ(x.data());
    Eigen::Matrix<double, 4, 3> J;
    J <<  0.5*xQ.w(), -0.5*xQ.z(),  0.5*xQ.y(),
          0.5*xQ.z(),  0.5*xQ.w(), -0.5*xQ.x(),
         -0.5*xQ.y(),  0.5*xQ.x(),  0.5*xQ.w(),
         -0.5*xQ.x(), -0.5*xQ.y(), -0.5*xQ.z();
    return J;
  }

  void ExpMapQuaternion::applyDiffRetractation_(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x, ReusableTemporaryMap& m)
  {
    pgs_assert(in.cols() == OutputDim_ && "Dimensions mismatch" );
    Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> a = m.getMap(in.rows(),3);
    a.noalias() = in*diffRetractation_(x);
    out = a;
  }

  Eigen::Matrix<double, 3, 4> ExpMapQuaternion::diffPseudoLog0_(const ConstRefVec& v)
  {
    const toConstQuat vQ(v.data());
    double n2 = vQ.vec().squaredNorm();
    double n = sqrt(n2);
    Eigen::Matrix<double, 3, 4> J;

    if (n < prec && vQ.w()!=0 )
    {
      double a = 2/vQ.w();
      double b = -2/(vQ.w()*vQ.w());
      J <<  a, 0, 0, b*vQ.x(),
            0, a, 0, b*vQ.y(),
            0, 0, a, b*vQ.z();
    }
    else
    {
      // log(x,y,z,w) = f(x,y,z,w)*[x;y;z]
      double f = atan2(2 * n * vQ.w(), vQ.w() * vQ.w() - n2) / n; 
      // df/dx = (x*atan((2*w*(x^2 + y^2 + z^2)^(1/2))/(- w^2 + x^2 + y^2 + z^2)))/(x^2 + y^2 + z^2)^(3/2) - ((2*w*x)/((x^2 + y^2 + z^2)^(1/2)*(- w^2 + x^2 + y^2 + z^2)) - (4*w*x*(x^2 + y^2 + z^2)^(1/2))/(- w^2 + x^2 + y^2 + z^2)^2)/(((4*w^2*(x^2 + y^2 + z^2))/(- w^2 + x^2 + y^2 + z^2)^2 + 1)*(x^2 + y^2 + z^2)^(1/2))
      // df/dy = (y*atan((2*w*(x^2 + y^2 + z^2)^(1/2))/(- w^2 + x^2 + y^2 + z^2)))/(x^2 + y^2 + z^2)^(3/2) - ((2*w*y)/((x^2 + y^2 + z^2)^(1/2)*(- w^2 + x^2 + y^2 + z^2)) - (4*w*y*(x^2 + y^2 + z^2)^(1/2))/(- w^2 + x^2 + y^2 + z^2)^2)/(((4*w^2*(x^2 + y^2 + z^2))/(- w^2 + x^2 + y^2 + z^2)^2 + 1)*(x^2 + y^2 + z^2)^(1/2))
      // df/dz = (z*atan((2*w*(x^2 + y^2 + z^2)^(1/2))/(- w^2 + x^2 + y^2 + z^2)))/(x^2 + y^2 + z^2)^(3/2) - ((2*w*z)/((x^2 + y^2 + z^2)^(1/2)*(- w^2 + x^2 + y^2 + z^2)) - (4*w*z*(x^2 + y^2 + z^2)^(1/2))/(- w^2 + x^2 + y^2 + z^2)^2)/(((4*w^2*(x^2 + y^2 + z^2))/(- w^2 + x^2 + y^2 + z^2)^2 + 1)*(x^2 + y^2 + z^2)^(1/2))
      
      // g = (atan((2*w*n)/(- w^2 + n2)))/(n2)^(3/2) - ((2*w)/(n*(- w^2 + n2)) - (4*w*n)/(- w^2 + n2)^2)/(((4*w^2*n2)/(- w^2 + n2)^2 + 1)*n)
      // g = (atan((2*w*n)/(- w^2 + n2)))/(n2*n) - ((2*w/n2)*(- w^2 + n2) - 4*w)/(4*w^2*n2+(- w^2 + n2)^2)
      // df/dx = g*x = (x*atan((2*w*n)/(- w^2 + n2)))/(n2)^(3/2) - ((2*w*x)/(n*(- w^2 + n2)) - (4*w*x*n)/(- w^2 + n2)^2)/(((4*w^2*n2)/(- w^2 + n2)^2 + 1)*n)
      // df/dy = g*y = (y*atan((2*w*n)/(- w^2 + n2)))/(n2)^(3/2) - ((2*w*y)/(n*(- w^2 + n2)) - (4*w*y*n)/(- w^2 + n2)^2)/(((4*w^2*n2)/(- w^2 + n2)^2 + 1)*n)
      // df/dz = g*z = (z*atan((2*w*n)/(- w^2 + n2)))/(n2)^(3/2) - ((2*w*z)/(n*(- w^2 + n2)) - (4*w*z*n)/(- w^2 + n2)^2)/(((4*w^2*n2)/(- w^2 + n2)^2 + 1)*n)
      // df/dw = -2/(w²+x²+y²+z²)
      double g = (atan((2*vQ.w()*n)/(- vQ.w()*vQ.w() + n2)))/(n2*n) - ((2*vQ.w()/n2)*(- vQ.w()*vQ.w() + n2) - 4*vQ.w())/(4*vQ.w()*vQ.w()*n2+(- vQ.w()*vQ.w() + n2)*(- vQ.w()*vQ.w() + n2));
      double dfdw = -2/(vQ.w()*vQ.w() + n2);
      /*
       * J = [ g.x²+f, g.y.x, g.z.x, df/dw.x]
       *     [ g.x.y, g.y²+f, g.z.y, df/dw.y]
       *     [ g.x.z, g.y.z, g.z²+f, df/dw.z]
      */
      J << g*vQ.x()*vQ.x()+f, g*vQ.y()*vQ.x(), g*vQ.z()*vQ.x(), dfdw*vQ.x(),
           g*vQ.x()*vQ.y(), g*vQ.y()*vQ.y()+f, g*vQ.z()*vQ.y(), dfdw*vQ.y(),
           g*vQ.x()*vQ.z(), g*vQ.y()*vQ.z(), g*vQ.z()*vQ.z()+f, dfdw*vQ.z();
    }
    return J;
  }

  void ExpMapQuaternion::applyDiffPseudoLog0_(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x, ReusableTemporaryMap& m)
  {
    pgs_assert(in.cols() == InputDim_ && "Dimensions mismatch" );
    Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> a = m.getMap(in.rows(),OutputDim_);
    a.noalias() = in*diffPseudoLog0_(x);
    out = a;
  }

  //void ExpMapQuaternion::applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec& v, ReusableTemporaryMap& m)
  void ExpMapQuaternion::applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec& , ReusableTemporaryMap& )
  {
    //OutputType E;
    //exponential(E,v);
    //Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> a = m.getMap(InputDim_, in.cols());
    //a.noalias() = E*in;
    //out = a;
    //TODO Make sure that out=in is OK here...
    out = in;
  }

  //void ExpMapQuaternion::applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec& v, ReusableTemporaryMap& m)
  void ExpMapQuaternion::applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec&, const ConstRefVec& , ReusableTemporaryMap&)
  {
    //OutputType E;
    //exponential(E,v);
    //Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> a = m.getMap(in.rows(), InputDim_);
    //a.noalias() = in*(E.transpose());
    //out = a;
    //TODO Make sure that out=in is OK here...
    out = in;
  }

  void ExpMapQuaternion::tangentConstraint_(RefMat, const ConstRefVec&)
  {
    //out is 0xt, no need to fill it
  }

  bool ExpMapQuaternion::isInTxM_(const ConstRefVec&, const ConstRefVec&, const double& )
  {
    return true;
  }

  void ExpMapQuaternion::forceOnTxM_(RefVec out, const ConstRefVec& in, const ConstRefVec&)
  {
    out = in;
  }

}


