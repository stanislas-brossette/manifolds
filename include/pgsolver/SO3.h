#ifndef _PGS_SO3_H_
#define _PGS_SO3_H_

#include <pgsolver/defs.h>
#include <pgsolver/Manifold.h>
#include <pgsolver/ReusableTemporaryMap.h>

namespace pgs
{
  template<typename Map>
  class SO3: public Manifold
  {
  public:
    SO3();
    virtual size_t numberOfSubmanifolds() const;
    virtual const Manifold& operator()(size_t i) const;
    virtual std::string toString(const ConstRefVec& val, const std::string& prefix = "") const;
  
  protected:
    //map operations
    virtual bool isValidInit_(const Eigen::VectorXd& val) const;
    virtual void plus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const;
    virtual void invMap_(RefVec out, const ConstRefVec& x) const;
    virtual void setIdentity_(RefVec out) const;
    virtual Eigen::MatrixXd diffMap_(const ConstRefVec& x) const;
    virtual void applyDiffMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    virtual Eigen::MatrixXd diffInvMap_(const ConstRefVec& x) const;
    virtual void applyDiffInvMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    virtual void applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const;

    mutable ReusableTemporaryMap bufferMap_;
  };

  //Implementations of the methods
  template<typename Map>
  inline SO3<Map>::SO3()
    : Manifold(Map::InputDim_, Map::OutputDim_)
  {
  }
  
  template<typename Map>
  inline bool SO3<Map>::isValidInit_(const Eigen::VectorXd& val) const
  {
    return Map::isValidInit_(val);
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
  inline std::string SO3<Map>::toString(const ConstRefVec& val, const std::string& prefix) const
  {
    std::string matPrefix = prefix + '[';
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", matPrefix, "]");
    std::stringstream ss;
    ss << (Eigen::Map<const typename Map::DisplayType>(val.data())).format(CleanFmt);
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
  inline void SO3<Map>::invMap_(RefVec out, const ConstRefVec& x) const
  {
    Map::invMap_( out, x);
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
    Map::applyDiffMap_(out, in, x, bufferMap_);
  }

  template<typename Map>
  inline Eigen::MatrixXd SO3<Map>::diffInvMap_(const ConstRefVec& x) const
  {
    return Map::diffInvMap_(x);
  }

  template<typename Map>
  inline void SO3<Map>::applyDiffInvMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    Map::applyDiffInvMap_(out, in, x, bufferMap_);
  }

  template<typename Map>
  inline void SO3<Map>::applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    Map::applyTransport_(out, in, x, v);
  }

  template<typename Map>
  inline void SO3<Map>::applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    Map::applyInvTransport_(out, in, x, v);
  }
}
#endif //_PGS_SO3_H_

