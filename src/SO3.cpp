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
    out = x + v;
  }

  void SO3::setIdentity_(Eigen::Ref<Eigen::VectorXd> out) const
  {
    out << 1,0,0,0,1,0,0,0,1;
  }
}
