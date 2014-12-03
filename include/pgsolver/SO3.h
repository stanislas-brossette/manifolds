#ifndef _PGS_SO3_H_
#define _PGS_SO3_H_

#include <pgsolver/defs.h>
#include <pgsolver/Manifold.h>

namespace pgs
{
  template<typename Map>
  class SO3: public Manifold
  {
  public:
    SO3();
    virtual bool isValidInit(const Eigen::VectorXd& val) const;
    virtual size_t numberOfSubmanifolds() const;
    virtual const Manifold& operator()(size_t i) const;
    virtual Segment getValue(RefVec val, size_t i) const;
    virtual ConstSegment getValueConst(const ConstRefVec& val, size_t i) const;
    virtual std::string toString(const ConstRefVec& val, std::string& prefix) const;
  protected:
    //map operations
    void plus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const;
    virtual void setIdentity_(RefVec out) const;
    virtual Eigen::MatrixXd diffMap_(const ConstRefVec& x) const;
    virtual void applyDiffMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
  };

  //Implementations of the methods
  template<typename Map>
  inline SO3<Map>::SO3()
    : Manifold(Map::InputDim_, Map::OutputDim_)
  {
  }
  
  template<typename Map>
  inline bool SO3<Map>::isValidInit(const Eigen::VectorXd& val) const
  {
    return Map::isValidInit(val);
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
  inline Segment SO3<Map>::getValue(RefVec val, size_t i) const
  {
    assert(i < 1 && "invalid index");
    return val.segment(0, static_cast<long> (representationDim()));
  }

  template<typename Map>
  inline ConstSegment SO3<Map>::getValueConst(const ConstRefVec& val, size_t i) const
  {
    assert(i < 1 && "invalid index");
    return val.segment(0,static_cast<long> (representationDim()));
  }
   
  template<typename Map>
  inline std::string SO3<Map>::toString(const ConstRefVec& val, std::string& prefix) const
  {
    std::string matPrefix = prefix + '[';
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", matPrefix, "]");
    std::stringstream ss;
    ss << (Eigen::Map<const typename Map::DisplayType>(val.data())).format(CleanFmt) << std::endl;
    return ss.str();
  }

  template<typename Map>
  inline void SO3<Map>::plus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    Map::plus_(out, x, v);
  }

  template<typename Map>
  inline void SO3<Map>::minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    Map::minus_( out, x, y);
  }
  
  template<typename Map>
  inline void SO3<Map>::setIdentity_(RefVec out) const
  {
    Map::setIdentity_(out);
  }

  template<typename Map>
  inline Eigen::MatrixXd SO3<Map>::diffMap_(const ConstRefVec& x) const
  {
    return Map::diffMap_(x);
  }

  template<typename Map>
  inline void SO3<Map>::applyDiffMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    Map::applyDiffMap_(out, in, x);
  }

}
#endif //_PGS_SO3_H_

