#ifndef _PGS_MANIFOLD_H_
#define _PGS_MANIFOLD_H_

#include "Point.h"

#include "defs.h"

namespace pgs
{
  class Manifold
  {
  public:
    Manifold(size_t dimension, size_t representationDimension);
    

    Point createPoint() const;
    Point createPoint(const Eigen::VectorXd& val) const;
    Point getIdentity() const;

    size_t dim() const;
    size_t representationDim() const;
    virtual size_t numberOfSubmanifolds() const = 0;
    virtual const Manifold& operator()(size_t i) const = 0;

    virtual Segment getValue(Eigen::Ref<Eigen::VectorXd> val, size_t i) const = 0;
    virtual ConstSegment getValueConst(const Eigen::Ref<const Eigen::VectorXd>& val, size_t i) const = 0;

    //map operations
    void setIdentity(Eigen::Ref<Eigen::VectorXd> out) const;
    void plus(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v) const;

    //for internal use
    void lock() const;

  protected:
    void setDimension(size_t d);
    void setRepresentationDimension(size_t rd);


    virtual void plus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v) const = 0;
    virtual void setIdentity_(Eigen::Ref<Eigen::VectorXd> out) const = 0;

    void testLock() const;

  private:
    size_t dimension_;
    size_t representationDim_;
    mutable bool lock_;
  };
}

#endif //_PGS_MANIFOLD_H_

