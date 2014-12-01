#include <cmath>
#include <pgsolver/SO3.h>

namespace pgs
{
  SO3::SO3()
    : Manifold(3, 9)
  {
  }

  size_t SO3::numberOfSubmanifolds() const
  {
    return 1;
  }
  const Manifold& SO3::operator()(size_t i) const
  {
    assert(i < 1 && "invalid index");
    return *this;
  }

  ConstSegment SO3::getValueConst(const Eigen::Ref<const Eigen::VectorXd>& val, size_t i) const
  {
    assert(i < 1 && "invalid index");
    return val.segment(0,static_cast<long> (representationDim()));
  }

  std::string SO3::toString(const Eigen::Ref<const Eigen::VectorXd>& val, std::string& prefix) const
  {
    std::string matPrefix = prefix + '[';
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", matPrefix, "]");
    std::stringstream ss;
    ss << (Eigen::Map<const Eigen::Matrix3d>(val.data())).format(CleanFmt) << std::endl;
    return ss.str();
  }

  Segment SO3::getValue(Eigen::Ref<Eigen::VectorXd> val, size_t i) const
  {
    assert(i < 1 && "invalid index");
    return val.segment(0, static_cast<long> (representationDim()));
  }

  void SO3::plus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v) const
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
    Eigen::Matrix3d E;
    E <<  1 - c*(v.y()*v.y() + v.z()*v.z()), -s*v.z() + c * v.x()*v.y(), s*v.y() + c * v.x()*v.z(),
          s * v.z() + c*v.x()*v.y(), 1 - c*(v.x()*v.x() + v.z()*v.z()) , -s*v.x() + c*v.y()*v.z(),
          -s*v.y() + c*v.x()*v.z(), s*v.x() + c*v.y()*v.z(), 1 - c*(v.x()*v.x() + v.y()*v.y());

    Eigen::Matrix3d rot;
    rot = (Eigen::Map<const Eigen::Matrix3d>(x.data()))*E;
    out = (Eigen::Map<const Eigen::VectorXd>(rot.data(),9));
  }
  
  void SO3::minus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y) const
  {
    Eigen::Matrix3d R(((Eigen::Map<const Eigen::Matrix3d>(y.data())).transpose())*(Eigen::Map<const Eigen::Matrix3d>(x.data())));
    Eigen::Vector3d v(-R(1,2), R(0,2), -R(0,1)); 
    double acosTr = std::acos((R.trace()-1)/2);
    if (v.norm() < 1e-8)
      out = v;
    else 
      {
        Eigen::Matrix3d diff(R-R.transpose());
        R = acosTr/(2*std::sin(acosTr))*(diff);
        v(0)=R(2,1);
        v(1)=R(0,2);
        v(2)=R(1,0);
        out = v;
      }
  }

  void SO3::setIdentity_(Eigen::Ref<Eigen::VectorXd> out) const
  {
    out << 1,0,0,0,1,0,0,0,1;
  }
}
