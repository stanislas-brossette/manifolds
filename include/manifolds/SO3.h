#ifndef _MANIFOLDS_SO3_H_
#define _MANIFOLDS_SO3_H_

#include <manifolds/defs.h>
#include <manifolds/Manifold.h>
#include <manifolds/ReusableTemporaryMap.h>

namespace pgs
{
  /// \brief Manifold representing the space of 3-dimensional rotations, also
  /// known as SO(3). It is templated by its map
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
    virtual bool isInM_(const Eigen::VectorXd& val) const;
    virtual void retractation_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void pseudoLog_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const;
    virtual void pseudoLog0_(RefVec out, const ConstRefVec& x) const;
    virtual void setZero_(RefVec out) const;
    virtual Eigen::MatrixXd diffRetractation_(const ConstRefVec& x) const;
    virtual void applyDiffRetractation_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    virtual Eigen::MatrixXd diffPseudoLog0_(const ConstRefVec& x) const;
    virtual void applyDiffPseudoLog0_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    virtual void applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const;

    virtual void tangentConstraint_(RefMat out, const ConstRefVec& x) const;
    virtual bool isInTxM_(const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void forceOnTxM_(RefVec out, const ConstRefVec& in, const ConstRefVec& x) const;
    virtual void limitMap_(RefVec out) const;

    mutable ReusableTemporaryMap bufferMap_;
  };

  //Implementations of the methods
  template<typename Map>
  inline SO3<Map>::SO3()
    : Manifold(3, Map::InputDim_, Map::OutputDim_)
  {
    name() = "SO3";
  }

  template<typename Map>
  inline bool SO3<Map>::isInM_(const Eigen::VectorXd& val) const
  {
    return Map::isInM_(val);
  }

  template<typename Map>
  inline size_t SO3<Map>::numberOfSubmanifolds() const
  {
    return 1;
  }

  template<typename Map>
  inline const Manifold& SO3<Map>::operator()(size_t i) const
  {
    pgs_assert(i < 1 && "invalid index");
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
  inline void SO3<Map>::retractation_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    Map::retractation_(out, x, v);
  }

  template<typename Map>
  inline void SO3<Map>::pseudoLog_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    Map::pseudoLog_( out, x, y);
  }

  template<typename Map>
  inline void SO3<Map>::pseudoLog0_(RefVec out, const ConstRefVec& x) const
  {
    Map::pseudoLog0_( out, x);
  }

  template<typename Map>
  inline void SO3<Map>::setZero_(RefVec out) const
  {
    Map::setZero_(out);
  }

  template<typename Map>
  inline Eigen::MatrixXd SO3<Map>::diffRetractation_(const ConstRefVec& x) const
  {
    return Map::diffRetractation_(x);
  }

  template<typename Map>
  inline void SO3<Map>::applyDiffRetractation_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    Map::applyDiffRetractation_(out, in, x, bufferMap_);
  }

  template<typename Map>
  inline Eigen::MatrixXd SO3<Map>::diffPseudoLog0_(const ConstRefVec& x) const
  {
    return Map::diffPseudoLog0_(x);
  }

  template<typename Map>
  inline void SO3<Map>::applyDiffPseudoLog0_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    Map::applyDiffPseudoLog0_(out, in, x, bufferMap_);
  }

  template<typename Map>
  inline void SO3<Map>::applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    Map::applyTransport_(out, in, x, v, bufferMap_);
  }

  template<typename Map>
  inline void SO3<Map>::applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    Map::applyInvTransport_(out, in, x, v, bufferMap_);
  }

  template<typename Map>
  void SO3<Map>::tangentConstraint_(RefMat out, const ConstRefVec& x) const
  {
    Map::tangentConstraint_(out, x);
  }

  template<typename Map>
  bool SO3<Map>::isInTxM_(const ConstRefVec& x, const ConstRefVec& v) const
  {
    return Map::isInTxM_(x, v);
  }

  template<typename Map>
  void SO3<Map>::forceOnTxM_(RefVec out, const ConstRefVec& in, const ConstRefVec& x) const
  {
    Map::forceOnTxM_(out, in, x);
  }

  template<typename Map>
  void SO3<Map>::limitMap_(RefVec out) const
  {
    double limitExpMap = M_PI/sqrt(3);
    out.setConstant(limitExpMap);
  }
}
#endif //_MANIFOLDS_SO3_H_

