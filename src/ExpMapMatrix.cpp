
#include <iostream>
#include <Eigen/Dense>
#include <pgsolver/defs.h>
#include <pgsolver/ExpMapMatrix.h>

namespace pgs
{
  const double ExpMapMatrix::prec = 1e-12;
 
  void ExpMapMatrix::plus_(RefVec out, ConstRefVec& x, ConstRefVec& v)
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

  void ExpMapMatrix::minus_(RefVec out, ConstRefVec& x, ConstRefVec& y)
  {
    typedef Eigen::Map<const Eigen::Matrix3d> ConstMapMat3;
    DisplayType R(((ConstMapMat3(y.data())).transpose())*(ConstMapMat3(x.data())));
    Eigen::Vector3d v(-R(1,2), R(0,2), -R(0,1)); 
    double acosTr = std::acos((R.trace()-1)/2);
    if (v.norm() < prec)
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
    out << 1,0,0,
           0,1,0,
           0,0,1;
  }

  bool ExpMapMatrix::isValidInit(const Eigen::VectorXd& val)
  {
    typedef Eigen::Map<const Eigen::Matrix3d> toMat3;
    bool out(val.size()==9);
    double det = toMat3(val.data()).determinant();
    out = out && (fabs(det - 1) < prec);
    out = out && 
      ((toMat3(val.data()).transpose())*toMat3(val.data())).isIdentity(prec);
    return out;
  }

  Eigen::MatrixXd ExpMapMatrix::diffMap_(ConstRefVec& )
  {
    Eigen::MatrixXd J(9,3);
    J << 0, 0, 0,  
         0, 0, 1,
         0,-1, 0,
         0, 0,-1,
         0, 0, 0,
         1, 0, 0,
         0, 1, 0,
        -1, 0, 0,
         0, 0, 0;
    return J;
  }

  void ExpMapMatrix::applyDiffMap_(RefMat , ConstRefVec& )
  {
  }
}

