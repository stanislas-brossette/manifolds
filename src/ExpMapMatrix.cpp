
#include <iostream>
#include <Eigen/Dense>
#include <pgsolver/defs.h>
#include <pgsolver/ExpMapMatrix.h>

namespace pgs
{
  int ExpMapMatrix::OutputDim() {return 9;}
  
  void ExpMapMatrix::plus_(RefVec out, ConstRefVec& x, ConstRefVec& v)
  {
    assert(v.size() == 3 && "Increment for expMap must be of size 3");
    double n = v.squaredNorm();
    assert(sqrt(n) < M_PI && "Increment for expMap must be of norm at most pi");
    double c, s;
    if (n < 1e-8)
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
                                                                                                                                                 
    DisplayType rot;
    rot = (Eigen::Map<const DisplayType>(x.data()))*E;
    out = (Eigen::Map<const OutputType>(rot.data(),9));
  }

  void ExpMapMatrix::minus_(RefVec out, ConstRefVec& x, ConstRefVec& y)
  {
    typedef Eigen::Map<const Eigen::Matrix3d> ConstMapMat3;
    DisplayType R(((ConstMapMat3(y.data())).transpose())*(ConstMapMat3(x.data())));
    Eigen::Vector3d v(-R(1,2), R(0,2), -R(0,1)); 
    double acosTr = std::acos((R.trace()-1)/2);
    if (v.norm() < 1e-8)
      out = v;
    else 
      {
        DisplayType diff(R-R.transpose());
        R = acosTr/(2*std::sin(acosTr))*(diff);
        v(0)=R(2,1);
        v(1)=R(0,2);
        v(2)=R(1,0);
        out = v;
      }
  }

  void ExpMapMatrix::setIdentity_(RefVec out)
  {
    out << 1,0,0,0,1,0,0,0,1;
  }

  bool ExpMapMatrix::isValidInit(const Eigen::VectorXd& val)
  {
    typedef Eigen::Map<const Eigen::Matrix3d> toMat3;
    bool out(val.size()==9);
    double det = toMat3(val.data()).determinant();
    Eigen::Matrix3d M = (((toMat3(val.data()).transpose())*toMat3(val.data()))
                              - Eigen::Matrix3d::Identity()).array().abs();
    out = out && (det - 1 < 1e-8);
    out = out && (M(0,0) <1.0e-8);
    out = out && (M(0,1) <1.0e-8);
    out = out && (M(0,2) <1.0e-8);
    out = out && (M(1,0) <1.0e-8);
    out = out && (M(1,1) <1.0e-8);
    out = out && (M(1,2) <1.0e-8);
    out = out && (M(2,0) <1.0e-8);
    out = out && (M(2,1) <1.0e-8);
    out = out && (M(2,2) <1.0e-8);
    return out;
  }
}

