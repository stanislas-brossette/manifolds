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

    virtual size_t numberOfSubmanifolds() const;
    virtual const Manifold& operator()(size_t i) const;

    //virtual Segment getValue(RefVec val, size_t i) const;
    //virtual ConstSegment getValueConst(const ConstRefVec& val, size_t i) const;
    //virtual Segment getValueTangent(RefVec val, size_t i) const;
    //virtual ConstSegment getValueTangentConst(const ConstRefVec& val, size_t i) const;
    virtual std::string toString(const ConstRefVec& val, const std::string& prefix = "") const;
  
  protected:
    //map operations
    virtual bool isValidInit_(const Eigen::VectorXd& ) const;
    virtual void plus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const;
    virtual void setIdentity_(RefVec out) const;
    virtual Eigen::MatrixXd diffMap_(const ConstRefVec& x) const;
    virtual void applyDiffMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    virtual Eigen::MatrixXd diffInvMap_(const ConstRefVec& x) const;
    virtual void applyDiffInvMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
  };
}

#endif //_PGS_REAL_SPACE_H_

