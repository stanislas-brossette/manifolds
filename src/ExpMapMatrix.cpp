#define _PGS_ASSERT_THROW_
#include <iostream>
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
  const double ExpMapMatrix::prec = 1e-12;
 
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

  Eigen::MatrixXd ExpMapMatrix::diffMap_(const ConstRefVec& x)
  {
    Eigen::MatrixXd J(9,3);
    J << 0   , -x(6), x(3),  
         0   , -x(7), x(4),
         0   , -x(8), x(5),
         x(6),  0   , -x(0),
         x(7),  0   , -x(1),
         x(8),  0   , -x(2),
        -x(3),  x(0), 0,
        -x(4),  x(1), 0,
        -x(5),  x(2), 0;
    return J;
  }

  void ExpMapMatrix::applyDiffMap_(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x)
  {
    //out.col(0) is used as a buffer to avoid aliasing in case where in and out
    //are the same memory array.

    std::string message("memory overlapping");
    pgs_assert(!utility::areOverlappingData(out.col(0), in.col(1)), message);
    pgs_assert(!utility::areOverlappingData(out.col(0), in.col(3)), message);
    pgs_assert(!utility::areOverlappingData(out.col(0), in.col(5)), message);
    pgs_assert(!utility::areOverlappingData(out.col(0), in.col(6)), message);
    pgs_assert(!utility::areOverlappingData(out.col(0), in.col(7)), message);
    pgs_assert(!utility::areOverlappingData(out.col(2), out.col(0)), message);
    pgs_assert(!utility::areOverlappingData(out.col(2), in.col(6)), message);
    pgs_assert(!utility::areOverlappingData(out.col(2), in.col(5)), message);
    pgs_assert(!utility::areOverlappingData(out.col(2), in.col(7)), message);
    pgs_assert(!utility::areOverlappingData(out.col(1), in.col(5)), message);
    pgs_assert(!utility::areOverlappingData(out.col(1), in.col(7)), message);

    applyDiffMapNoAssert_(out, in, x);
  }

  void ExpMapMatrix::applyDiffMapNoAssert_(
      RefMat out, const ConstRefMat& in, const ConstRefVec&)
  {

    out.col(0) = in.col(2);
    out.col(2) =  in.col(1) - in.col(3);
    out.col(1) = -out.col(0) + in.col(6);
    out.col(0) =  in.col(5) - in.col(7);
  }
  Eigen::MatrixXd ExpMapMatrix::diffInvMap_(const ConstRefVec& R)
  {
    Eigen::MatrixXd J(3,9);
    J.setZero();
    Eigen::Vector3d v;
    v << -R(7), R(6), -R(3);  //Valid approximation of the log when v<<1
    double trR = R(0)+R(4)+R(8); //Trace of R;

    if (v.norm()<prec) //Probably should not use prec here
    {
      std::cout << "small angle: "<< std::endl;
      double f = (acos((trR-1)/2/M_PI))/sin(acos((trR-1)/2/M_PI));
      std::cout << "f = " << f << std::endl;
      double df = -1/6+1/15*(trR-3);
      std::cout << "df = " << df << std::endl;
      J.col(0) = df*v;
      J.col(1) = f*Eigen::Vector3d(0   , 0   , 1/2 );
      J.col(2) = f*Eigen::Vector3d(0   , -1/2, 0   );
      J.col(3) = f*Eigen::Vector3d(0   , 0   , -1/2);
      J.col(4) = df*v;
      J.col(5) = f*Eigen::Vector3d(1/2 , 0   , 0   );
      J.col(6) = f*Eigen::Vector3d(0   , 1/2 , 0   );
      J.col(7) = f*Eigen::Vector3d(-1/2, 0   , 0   );
      J.col(8) = df*v;
    }
    else
    {
      std::cout << "small angle: "<< std::endl;
      double trRm1o2 = (trR-1)/2;
      std::cout << "trRm1o2 = " << trRm1o2 << std::endl;
      double f = acos(trRm1o2)/sin(acos(trRm1o2));
      std::cout << "f = " << f << std::endl;
      double df = 1/(2*(pow(trRm1o2,2)-1))+(trRm1o2*acos(trRm1o2))/(2*pow((1-pow(trRm1o2,2)),(1.5)));
      std::cout << "df = " << df << std::endl;
      Eigen::Vector3d g(R(5) - R(7), R(6) - R(2), R(1) - R(3));
      std::cout << "g = " << g << std::endl;
      J.col(0) << 0.5*df*g;
      J.col(1) << 0 , 0 , 0.5*f ;
      J.col(2) << 0 , -0.5*f, 0 ;
      J.col(3) << 0 , 0 , -0.5*f;
      J.col(4) = J.col(0); 
      J.col(5) << 0.5*f , 0 , 0 ;
      J.col(6) << 0 , 0.5*f , 0 ;
      J.col(7) << -0.5*f, 0 , 0 ;
      J.col(8) = J.col(0); 
    }

    return J;
  }
}

