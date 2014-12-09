#ifndef _PGS_MANIFOLD_H_
#define _PGS_MANIFOLD_H_

#include <iostream>
#include <Eigen/Core>
#include <pgsolver/Point.h>
#include <pgsolver/defs.h>
#include <pgsolver/RefCounter.h>

namespace pgs
{
  enum eDimension
  {
    R,  //representation space
    T,  //tangent space
    F,  //full space
  };


  class Manifold : public RefCounter
  {
  public:
    Manifold(Index dimension, Index representationDimension);

    Point createPoint() const;
    Point createPoint(const Eigen::VectorXd& val) const;
    Point getIdentity() const;

    virtual bool isValidInit(const Eigen::VectorXd& val) const;

    Index dim() const;
    Index representationDim() const;
    virtual size_t numberOfSubmanifolds() const = 0;
    virtual const Manifold& operator()(size_t i) const = 0;

    //view
    Segment getValue(RefVec val, size_t i) const { return getView<R>(val, i); }
    ConstSegment getValueConst(const ConstRefVec& val, size_t i) const { return getView<R>(val, i); }
    Segment getValueTangent(RefVec val, size_t i) const { return getView<T>(val, i); }
    ConstSegment getValueTangentConst(const ConstRefVec& val, size_t i) const { return getView<T>(val, i); }

    template<int D> Segment getView(RefVec val, size_t i) const;
    template<int D> ConstSegment getView(const ConstRefVec& val, size_t i) const;

    template<int Dr, int Dc> Eigen::MatrixXd getView(RefMat val, size_t i) const;
    template<int Dr, int Dc> Eigen::MatrixXd getView(const ConstRefMat& val, size_t i) const;
    template<int Dr, int Dc> Eigen::MatrixXd getView(RefMat val, size_t ir, size_t ic) const;
    template<int Dr, int Dc> Eigen::MatrixXd getView(const ConstRefMat& val, size_t ir, size_t ic) const;

    virtual std::string toString(const ConstRefVec& val, std::string& prefix) const = 0;

    //map operations
    void setIdentity(RefVec out) const;
    void plus(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const;
    void minus(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const;
    void invMap(RefVec out, const ConstRefVec& x) const;
    Eigen::MatrixXd diffMap(const ConstRefVec& x) const;
    void applyDiffMap(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    Eigen::MatrixXd diffInvMap(const ConstRefVec& x) const;
    void applyDiffInvMap(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;

    //for internal use
    void lock() const;

  protected:
    void setDimension(Index d);
    void setRepresentationDimension(Index rd);
    virtual bool isValidInit_(const Eigen::VectorXd& val) const = 0;

    template<int D>
    Index getDim() const;
    template<int D>
    Index getDim(size_t i) const;
    template<int D>
    Index getStart(size_t i) const;
    virtual Index startR(size_t i) const;
    virtual Index startT(size_t i) const;

    virtual void plus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const = 0;
    virtual void minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const = 0;
    virtual void setIdentity_(RefVec out) const = 0;
    virtual Eigen::MatrixXd diffMap_(const ConstRefVec& x) const = 0;
    virtual void applyDiffMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const = 0;
    virtual Eigen::MatrixXd diffInvMap_(const ConstRefVec& x) const = 0;
    virtual void applyDiffInvMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const = 0;

    void testLock() const;

  private:
    Index dimension_;
    Index representationDim_;
    mutable bool lock_;
  };

  template<int D>
  inline Segment Manifold::getView(RefVec val, size_t i) const
  {
    assert(i < numberOfSubmanifolds() && "invalid index");
    assert(val.size() == getDim<D>());
    return val.segment(getStart<D>(i), getDim<D>(i));
  }
  
  template<int D>
  inline ConstSegment Manifold::getView(const ConstRefVec& val, size_t i) const
  {
    assert(i < numberOfSubmanifolds() && "invalid index");
    assert(val.size() == getDim<D>());
    return val.segment(getStart<D>(i), getDim<D>(i));
  }

  template<>
  inline Index Manifold::getStart<R>(size_t i) const
  {
    assert(i < numberOfSubmanifolds() && "invalid index");
    return startR(i);
  }

  template<>
  inline Index Manifold::getStart<T>(size_t i) const
  {
    assert(i < numberOfSubmanifolds() && "invalid index");
    return startT(i);
  }

  template<>
  inline Index Manifold::getDim<R>() const
  {
    return representationDim();
  }

  template<>
  inline Index Manifold::getDim<T>() const
  {
    return dim();
  }
  
  template<>
  inline Index Manifold::getDim<R>(size_t i) const
  {
    assert(i < numberOfSubmanifolds() && "invalid index");
    return operator()(i).representationDim();
  }

  template<>
  inline Index Manifold::getDim<T>(size_t i) const
  {
    assert(i < numberOfSubmanifolds() && "invalid index");
    return operator()(i).dim();
  }

  inline Index Manifold::startR(size_t i) const
  {
    assert(i < numberOfSubmanifolds() && "invalid index");
    return 0;
  }

  inline Index Manifold::startT(size_t i) const
  {
    assert(i < numberOfSubmanifolds() && "invalid index");
    return 0;
  }


  template<int Dr, int Dc>
  inline Eigen::MatrixXd Manifold::getView(RefMat val, size_t i) const
  {
    return getView<Dr, Dc>(val, i, i);
  }

  template<int Dr, int Dc>
  inline Eigen::MatrixXd Manifold::getView(const ConstRefMat& val, size_t i) const
  {
    return getView<Dr, Dc>(val, i, i);
  }

  template<int Dr, int Dc>
  inline Eigen::MatrixXd Manifold::getView(RefMat val, size_t ir, size_t ic) const
  {
    return val.block(getStart<Dr>(ir),
      getStart<Dc>(ic),
      getDim<Dr>(ir),
      getDim<Dc>(ir));
  }

  template<int Dr, int Dc>
  inline Eigen::MatrixXd Manifold::getView(const ConstRefMat& val, size_t ir, size_t ic) const
  {
    return val.block(getStart<Dr>(ir),
      getStart<Dc>(ic),
      getDim<Dr>(ir),
      getDim<Dc>(ir));
  }

  template<>
  inline Eigen::MatrixXd Manifold::getView<F, R>(RefMat val, size_t, size_t ic) const
  {
    return val.middleCols(getStart<R>(ic), getDim<R>(ic));
  }

  template<>
  inline Eigen::MatrixXd Manifold::getView<F, R>(const ConstRefMat& val, size_t, size_t ic) const
  {
    return val.middleCols(getStart<R>(ic), getDim<R>(ic));
  }

  template<>
  inline Eigen::MatrixXd Manifold::getView<F, T>(RefMat val, size_t, size_t ic) const
  {
    return val.middleCols(getStart<T>(ic), getDim<T>(ic));
  }

  template<>
  inline Eigen::MatrixXd Manifold::getView<F, T>(const ConstRefMat& val, size_t, size_t ic) const
  {
    return val.middleCols(getStart<T>(ic), getDim<T>(ic));
  }

  template<>
  inline Eigen::MatrixXd Manifold::getView<R, F>(RefMat val, size_t ir, size_t) const
  {
    return val.middleRows(getStart<R>(ir), getDim<R>(ir));
  }

  template<>
  inline Eigen::MatrixXd Manifold::getView<R, F>(const ConstRefMat& val, size_t ir, size_t) const
  {
    return val.middleRows(getStart<R>(ir), getDim<R>(ir));
  }

  template<>
  inline Eigen::MatrixXd Manifold::getView<T, F>(RefMat val, size_t ir, size_t) const
  {
    return val.middleRows(getStart<T>(ir), getDim<T>(ir));
  }

  template<>
  inline Eigen::MatrixXd Manifold::getView<T, F>(const ConstRefMat& val, size_t ir, size_t) const
  {
    return val.middleRows(getStart<T>(ir), getDim<T>(ir));
  }

  template<>
  inline Eigen::MatrixXd Manifold::getView<F, F>(RefMat val, size_t, size_t) const
  {
    return val;
  }

  template<>
  inline Eigen::MatrixXd Manifold::getView<F, F>(const ConstRefMat& val, size_t, size_t) const
  {
    return val;
  }

}

#endif //_PGS_MANIFOLD_H_

