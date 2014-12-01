#ifndef _PGS_CARTESIAN_PRODUCT_H_
#define _PGS_CARTESIAN_PRODUCT_H_

#include <vector>
#include <pgsolver/Manifold.h>

namespace pgs
{
  class CartesianProduct : public Manifold
  {
  public:
    CartesianProduct();
    CartesianProduct(const Manifold& m1, const Manifold& m2);

    CartesianProduct& multiply(const Manifold& m);

    virtual size_t numberOfSubmanifolds() const;
    virtual const Manifold& operator()(size_t i) const;

    virtual Segment getValue(Eigen::Ref<Eigen::VectorXd> val, size_t i) const;
    virtual ConstSegment getValueConst(const Eigen::Ref<const Eigen::VectorXd>& val, size_t i) const;
    virtual std::string toString(const Eigen::Ref<const Eigen::VectorXd>& val, std::string& prefix) const;
  


  protected:
    virtual void plus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& v) const;
    virtual void minus_(Eigen::Ref<Eigen::VectorXd> out, const Eigen::Ref<const Eigen::VectorXd>& x, const Eigen::Ref<const Eigen::VectorXd>& y) const;
    virtual void setIdentity_(Eigen::Ref<Eigen::VectorXd> out) const;

  private:
    std::vector<const Manifold* > submanifolds_;
    std::vector<Index> startIndexT_;
    std::vector<Index> startIndexR_;
  };
}

#endif //_PGS_CARTESIAN_PRODUCT_H_

