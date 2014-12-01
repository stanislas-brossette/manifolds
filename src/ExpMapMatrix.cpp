
#include <iostream>
#include <pgsolver/ExpMapMatrix.h>

namespace pgs
{
  int ExpMapMatrix::OutputDim() {return 9;}
  
  void ExpMapMatrix::plus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v)
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
    E <<  1 - c*(v.y()*v.y() + v.z()*v.z()), -s*v.z() + c * v.x()*v.y(), s*v.y() + c * v.x()*v.z(),
          s * v.z() + c*v.x()*v.y(), 1 - c*(v.x()*v.x() + v.z()*v.z()) , -s*v.x() + c*v.y()*v.z(),
          -s*v.y() + c*v.x()*v.z(), s*v.x() + c*v.y()*v.z(), 1 - c*(v.x()*v.x() + v.y()*v.y());
                                                                                                                                                 
    DisplayType rot;
    rot = (Eigen::Map<const DisplayType>(x.data()))*E;
    out = (Eigen::Map<const OutputType>(rot.data(),9));
  }

  void ExpMapMatrix::minus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y)
  {
    DisplayType R(((Eigen::Map<const Eigen::Matrix3d>(y.data())).transpose())*(Eigen::Map<const Eigen::Matrix3d>(x.data())));
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

  void ExpMapMatrix::setIdentity_(Eigen::Ref<Eigen::VectorXd> out)
  {
    out << 1,0,0,0,1,0,0,0,1;
  }

  bool ExpMapMatrix::isValidInit(const Eigen::VectorXd& val)
  {
    std::cout << "checkValidRotation called" << std::endl;
    std::cout << val.transpose() << std::endl;
    return true;
  }
}

