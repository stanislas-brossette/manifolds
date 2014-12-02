#ifndef _PGS_REAL_SPACE_H_
#define _PGS_REAL_SPACE_H_

#include <iostream>
#include <pgsolver/Manifold.h>

namespace pgs
{
  class RealSpace: public Manifold
  {
  public:
    RealSpace(Index n);

    virtual bool isValidInit(const Eigen::VectorXd& ) const;

    virtual size_t numberOfSubmanifolds() const;
    virtual const Manifold& operator()(size_t i) const;

    virtual Segment getValue(RefVec val, size_t i) const;
    virtual ConstSegment getValueConst(ConstRefVec& val, size_t i) const;
    virtual std::string toString(ConstRefVec& val, std::string& prefix) const;
  
  protected:
    //map operations
    virtual void plus_(RefVec out, ConstRefVec& x, ConstRefVec& v) const;
    virtual void minus_(RefVec out, ConstRefVec& x, ConstRefVec& y) const;

    virtual void setIdentity_(RefVec out) const;
  };
}

#endif //_PGS_REAL_SPACE_H_

