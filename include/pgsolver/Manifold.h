#ifndef _PGS_MANIFOLD_H_
#define _PGS_MANIFOLD_H_

#include <iostream>
#include <Eigen/Core>
#include <pgsolver/Point.h>
#include <pgsolver/defs.h>
#include <pgsolver/RefCounter.h>

namespace pgs
{
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

    virtual Segment getValue(RefVec val, size_t i) const = 0;
    virtual ConstSegment getValueConst(const ConstRefVec& val, size_t i) const = 0;
    virtual Segment getValueTangent(RefVec val, size_t i) const = 0;
    virtual ConstSegment getValueTangentConst(const ConstRefVec& val, size_t i) const = 0;

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
}

#endif //_PGS_MANIFOLD_H_

