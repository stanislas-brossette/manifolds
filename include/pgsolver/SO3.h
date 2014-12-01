#ifndef _PGS_SO3_H_
#define _PGS_SO3_H_

#include <pgsolver/Manifold.h>

namespace pgs
{
  template<typename Map>
  class SO3: public Manifold
  {
  public:
    SO3();
    virtual size_t numberOfSubmanifolds() const;
    virtual const Manifold& operator()(size_t i) const;
    virtual Segment getValue(Eigen::Ref<Eigen::VectorXd> val, size_t i) const;
    virtual ConstSegment getValueConst(const Eigen::Ref<const Eigen::VectorXd>& val, size_t i) const;
    virtual std::string toString(const Eigen::Ref<const Eigen::VectorXd>& val, std::string& prefix) const;
  protected:
    //map operations
    void plus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v) const;
    virtual void minus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y) const;
    virtual void setIdentity_(Eigen::Ref<Eigen::VectorXd> out) const;
  };


  template<typename Map>
  inline SO3<Map>::SO3()
    : Manifold(3, Map::OutputDim())
  {
  }

  template<typename Map>
  inline size_t SO3<Map>::numberOfSubmanifolds() const
  {
    return 1;
  }

  template<typename Map>
  inline const Manifold& SO3<Map>::operator()(size_t i) const
  {
    assert(i < 1 && "invalid index");
    return *this;
  }

  template<typename Map>
  inline Segment SO3<Map>::getValue(Eigen::Ref<Eigen::VectorXd> val, size_t i) const
  {
    assert(i < 1 && "invalid index");
    return val.segment(0, static_cast<long> (representationDim()));
  }

  template<typename Map>
  inline ConstSegment SO3<Map>::getValueConst(const Eigen::Ref<const Eigen::VectorXd>& val, size_t i) const
  {
    assert(i < 1 && "invalid index");
    return val.segment(0,static_cast<long> (representationDim()));
  }
   
  template<typename Map>
  inline std::string SO3<Map>::toString(const Eigen::Ref<const Eigen::VectorXd>& val, std::string& prefix) const
  {
    std::string matPrefix = prefix + '[';
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", matPrefix, "]");
    std::stringstream ss;
    ss << (Eigen::Map<const typename Map::DisplayType>(val.data())).format(CleanFmt) << std::endl;
    return ss.str();
  }

  template<typename Map>
  inline void SO3<Map>::plus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v) const
  {
    Map::plus_(out, x, v);
  }

  template<typename Map>
  inline void SO3<Map>::minus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y) const
  {
    Map::minus_( out, x, y);
  }
  
  template<typename Map>
  inline void SO3<Map>::setIdentity_(Eigen::Ref<Eigen::VectorXd> out) const
  {
    Map::setIdentity_(out);
  }
}
#endif //_PGS_SO3_H_

