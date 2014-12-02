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

    virtual bool isValidInit(const Eigen::VectorXd& val) const;

    CartesianProduct& multiply(const Manifold& m);

    virtual size_t numberOfSubmanifolds() const;
    virtual const Manifold& operator()(size_t i) const;

    virtual Segment getValue(RefVec val, size_t i) const;
    virtual ConstSegment getValueConst(ConstRefVec& val, size_t i) const;
    virtual std::string toString(ConstRefVec& val, std::string& prefix) const;
  


  protected:
    virtual void plus_(RefVec out, ConstRefVec& x, ConstRefVec& v) const;
    virtual void minus_(RefVec out, ConstRefVec& x, ConstRefVec& y) const;
    virtual void setIdentity_(RefVec out) const;
    virtual Eigen::MatrixXd diffMap_(ConstRefVec& x) const;
    virtual void applyDiffMap_(RefMat inOut, ConstRefVec& x) const;

  private:
    std::vector<const Manifold* > submanifolds_;
    std::vector<Index> startIndexT_;
    std::vector<Index> startIndexR_;
  };
}

#endif //_PGS_CARTESIAN_PRODUCT_H_

