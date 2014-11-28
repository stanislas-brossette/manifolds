#ifndef _PGS_MANIFOLD_H_
#define _PGS_MANIFOLD_H_

#include <iostream>
#include "Point.h"
#include "defs.h"

namespace pgs
{
  class Manifold
  {
  public:
    Manifold(Index dimension, Index representationDimension);
    

    Point createPoint() const;
    Point createPoint(const Eigen::VectorXd& val) const;
    Point getIdentity() const;

    Index dim() const;
    Index representationDim() const;
    virtual size_t numberOfSubmanifolds() const = 0;
    virtual const Manifold& operator()(size_t i) const = 0;

    virtual Segment getValue(Eigen::Ref<Eigen::VectorXd> val, size_t i) const = 0;
    virtual ConstSegment getValueConst(const Eigen::Ref<const Eigen::VectorXd>& val, size_t i) const = 0;

    virtual std::string toString(const Eigen::Ref<const Eigen::VectorXd>& val, std::string& prefix) const = 0;
    //map operations
    void setIdentity(Eigen::Ref<Eigen::VectorXd> out) const;
    void plus(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v) const;

    //for internal use
    void lock() const;

  protected:
    void setDimension(Index d);
    void setRepresentationDimension(Index rd);


    virtual void plus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v) const = 0;
    virtual void setIdentity_(Eigen::Ref<Eigen::VectorXd> out) const = 0;

    void testLock() const;

  private:
    Index dimension_;
    Index representationDim_;
    mutable bool lock_;
  };
}

#endif //_PGS_MANIFOLD_H_

