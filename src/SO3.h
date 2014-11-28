#ifndef _PGS_SO3_H_
#define _PGS_SO3_H_

#include "Manifold.h"

namespace pgs
{
  class SO3: public Manifold
  {
  public:
    SO3();

    virtual size_t numberOfSubmanifolds() const;
    virtual const Manifold& operator()(size_t i) const;

    virtual Segment getValue(Eigen::Ref<Eigen::VectorXd> val, size_t i) const;
    virtual ConstSegment getValueConst(const Eigen::Ref<const Eigen::VectorXd>& val, size_t i) const;
    virtual std::string toString(const Eigen::Ref<const Eigen::VectorXd>& val, std::string& prefix) const;
  


  protected:
    //map operations
    virtual void plus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v) const;

    virtual void setIdentity_(Eigen::Ref<Eigen::VectorXd> out) const;
  };
}

#endif //_PGS_SO3_H_

