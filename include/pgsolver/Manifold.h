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

    virtual bool isValidInit(const Eigen::VectorXd& val) const = 0;

    Index dim() const;
    Index representationDim() const;
    virtual size_t numberOfSubmanifolds() const = 0;
    virtual const Manifold& operator()(size_t i) const = 0;

    virtual Segment getValue(RefVec val, size_t i) const = 0;
    virtual ConstSegment getValueConst(ConstRefVec& val, size_t i) const = 0;

    virtual std::string toString(ConstRefVec& val, std::string& prefix) const = 0;
    //map operations
    void setIdentity(RefVec out) const;
    void plus(RefVec out, ConstRefVec& x, ConstRefVec& v) const;
    void minus(RefVec out, ConstRefVec& x, ConstRefVec& v) const;
    

    //for internal use
    void lock() const;

  protected:
    void setDimension(Index d);
    void setRepresentationDimension(Index rd);


    virtual void plus_(RefVec out, ConstRefVec& x, ConstRefVec& v) const = 0;
    virtual void minus_(RefVec out, ConstRefVec& x, ConstRefVec& v) const = 0;
    virtual void setIdentity_(RefVec out) const = 0;

    void testLock() const;

  private:
    Index dimension_;
    Index representationDim_;
    mutable bool lock_;
  };
}

#endif //_PGS_MANIFOLD_H_

