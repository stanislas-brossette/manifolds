#include <iostream>
#include <boost/math/special_functions/sinc.hpp>
#include <Eigen/Dense>
#include <pgsolver/defs.h>
#include <pgsolver/ExpMapMatrix.h>
#include <pgsolver/pgs_assert.h>

namespace utility
{
  // areOverlappingData tests if the data pointed by a and b are overlapping, in
  // which case, some aliasing could appear.
  // TODO, explain more.
  // This function assumes that there wasn't any copies of a or b before the
  // call.
  bool areOverlappingData(const pgs::ConstRefVec& a, const pgs::ConstRefVec& b)
  {
    bool res = false;
    if (&a.coeff(a.rows()-1) < &b.coeff(0) || &b.coeff(b.rows()-1) < &a.coeff(0))
    {
      return false;
    }

    for (int i = 0;i<a.rows();++i)
    {
      for (int j = i;j<b.rows();++j)
      {
        res = res || (&a.coeff(i) == &b.coeff(j));
      }
    }
    return res;
  }
}

namespace pgs
{
  const double ExpMapMatrix::prec = 1e-8;
 
  void ExpMapMatrix::plus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v)
  {
    assert(v.size() == 3 && "Increment for expMap must be of size 3");
    double n = v.squaredNorm();
    assert(sqrt(n) < M_PI && "Increment for expMap must be of norm at most pi");
    double c, s;
    if (n < prec)
    {
      c = 0.5 - n / 24;
      s = 1 - n / 6;
    }
    else
    {
      double t = sqrt(n);
      c = (1 - cos(t)) / n;
      s = sin(t) / t;
    }
    DisplayType E;
    E <<  1 - c*(v.y()*v.y() + v.z()*v.z()), 
          -s*v.z() + c * v.x()*v.y(), 
          s*v.y() + c * v.x()*v.z(),
          s * v.z() + c*v.x()*v.y(), 
          1 - c*(v.x()*v.x() + v.z()*v.z()) , 
          -s*v.x() + c*v.y()*v.z(),
          -s*v.y() + c*v.x()*v.z(), 
          s*v.x() + c*v.y()*v.z(), 
          1 - c*(v.x()*v.x() + v.y()*v.y());
                                                                                                                                                 
    Eigen::Map<DisplayType>(out.data()) = (Eigen::Map<const DisplayType>(x.data()))*E;
  }

  void ExpMapMatrix::minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& y)
  {
    typedef Eigen::Map<const Eigen::Matrix3d> ConstMapMat3;
    DisplayType R(((ConstMapMat3(y.data())).transpose())*(ConstMapMat3(x.data())));
    logarithm(out,R);
  }

  void ExpMapMatrix::invMap_(RefVec out, const ConstRefVec& x)
  {
    typedef Eigen::Map<const Eigen::Matrix3d> ConstMapMat3;
    DisplayType R(ConstMapMat3(x.data()));
    logarithm(out,R);
  }
  
  void ExpMapMatrix::logarithm(RefVec out, const DisplayType& R)
  {
    Eigen::Vector3d v(-R(1,2), R(0,2), -R(0,1)); 
    double acosTr = std::acos((R.trace()-1)/2);
    if (v.norm() < prec)
      out = v;
    else 
    {
      DisplayType diff(R-R.transpose());
      double coeff = acosTr/(2*std::sin(acosTr));
      v(0)=diff(2,1)*coeff;
      v(1)=diff(0,2)*coeff;
      v(2)=diff(1,0)*coeff;
      out = v;
    }
  }

  void ExpMapMatrix::setIdentity_(RefVec out)
  {
    out << 1,0,0,
           0,1,0,
           0,0,1;
  }

  bool ExpMapMatrix::isValidInit_(const Eigen::VectorXd& val)
  {
    typedef Eigen::Map<const Eigen::Matrix3d> toMat3;
    bool out(val.size()==9);
    double det = toMat3(val.data()).determinant();
    out = out && (fabs(det - 1) < prec);
    out = out && 
      ((toMat3(val.data()).transpose())*toMat3(val.data())).isIdentity(prec);
    return out;
  }

  Eigen::Matrix<double, 9, 3> ExpMapMatrix::diffMap_(const ConstRefVec& x)
  {
    Eigen::Matrix<double, 9, 3> J;
    J << 0   , -x(6), x(3) ,  
         0   , -x(7), x(4) ,
         0   , -x(8), x(5) ,
         x(6),  0   , -x(0),
         x(7),  0   , -x(1),
         x(8),  0   , -x(2),
        -x(3),  x(0), 0    ,
        -x(4),  x(1), 0    ,
        -x(5),  x(2), 0    ;
    return J;
  }

  void ExpMapMatrix::applyDiffMap_(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x, ReusableTemporaryMap& m)
  {
    assert(in.cols() == OutputDim_ && "Dimensions mismatch" );
    Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> a = m.getMap(in.rows(),3);
    a.noalias() = in*diffMap_(x);
    out = a;
  }

  Eigen::Matrix<double, 3, 9> ExpMapMatrix::diffInvMap_(const ConstRefVec& R)
  {
    Eigen::Matrix<double, 3, 9> J;
    J.setZero();
    //Valid approximation of the log when v<<1
    Eigen::Vector3d v((R(5)-R(7))/2, (R(6)-R(2))/2, (R(1)-R(3))/2);      

    double trR = R(0)+R(4)+R(8); //Trace of R;
    double trRm1o2 = (trR-1) / 2;
    double s2 = 1 - trRm1o2*trRm1o2;
    double f = 1/boost::math::sinc_pi(acos(trRm1o2)); //=acos(trRm1o2)/sin(acos(trRm1o2));
    double hf = f/2;
    double df;
    if (v.squaredNorm()<1e-10) //Probably should not use this switch value here
    {
      //TODO: I do not know how to test that. Verification with finite
      //differences proved unefficient. Need to get back to it later.
      //
      //Here we use the Taylor approximation of the diff of the log
      double x = trR-3;
      df = trR/15 - (3*(x*x))/140 + (2*(x*x*x))/315 - (5*(x*x*x*x))/2772 - 11/30;
    }
    else
    {
      df = (trRm1o2*f-1) / (2 * s2);
    }
    J.col(0) = df * v;
    J.col(1) = Eigen::Vector3d( 0 , 0 , hf);
    J.col(2) = Eigen::Vector3d( 0 ,-hf, 0 );
    J.col(3) = Eigen::Vector3d( 0 , 0 ,-hf);
    J.col(4) = J.col(0);
    J.col(5) = Eigen::Vector3d( hf, 0 , 0 );
    J.col(6) = Eigen::Vector3d( 0 , hf, 0 );
    J.col(7) = Eigen::Vector3d(-hf, 0 , 0 );
    J.col(8) = J.col(0);
    return J;
  }

  void ExpMapMatrix::applyDiffInvMap_(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x, ReusableTemporaryMap& m)
  {
    assert(in.cols() == InputDim_ && "Dimensions mismatch" );
    Eigen::Map<Eigen::MatrixXd, Eigen::Aligned> a = m.getMap(in.rows(),9);
    a.noalias() = in*diffInvMap_(x);
    out = a;
  }
}

