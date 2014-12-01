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

    virtual bool isValidInit(const Eigen::VectorXd& val) const;

    virtual size_t numberOfSubmanifolds() const;
    virtual const Manifold& operator()(size_t i) const;

    virtual Segment getValue(Eigen::Ref<Eigen::VectorXd> val, size_t i) const;
    virtual ConstSegment getValueConst(const Eigen::Ref<const Eigen::VectorXd>& val, size_t i) const;
    virtual std::string toString(const Eigen::Ref<const Eigen::VectorXd>& val, std::string& prefix) const;
  
  protected:
    //map operations
    virtual void plus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v) const;
    virtual void minus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y) const;

    virtual void setIdentity_(Eigen::Ref<Eigen::VectorXd> out) const;
  };
}

#endif //_PGS_REAL_SPACE_H_

