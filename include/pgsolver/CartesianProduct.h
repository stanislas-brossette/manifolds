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

    //virtual Segment getValue(RefVec val, size_t i) const;
    //virtual ConstSegment getValueConst(const ConstRefVec& val, size_t i) const;
    //virtual Segment getValueTangent(RefVec val, size_t i) const;
    //virtual ConstSegment getValueTangentConst(const ConstRefVec& val, size_t i) const;
    virtual std::string toString(const ConstRefVec& val, const std::string& prefix = "") const;
  


  protected:
    virtual bool isValidInit_(const Eigen::VectorXd& val) const;

    virtual Index startR(size_t i) const;
    virtual Index startT(size_t i) const;

    virtual void plus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const;
    virtual void setIdentity_(RefVec out) const;
    virtual Eigen::MatrixXd diffMap_(const ConstRefVec& x) const;
    virtual void applyDiffMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    virtual Eigen::MatrixXd diffInvMap_(const ConstRefVec& x) const;
    virtual void applyDiffInvMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;

  private:
    std::vector<const Manifold* > submanifolds_;
    std::vector<Index> startIndexT_;
    std::vector<Index> startIndexR_;
  };

  inline Index CartesianProduct::startR(size_t i) const
  {
    assert(i < numberOfSubmanifolds() && "invalid index");
    return startIndexR_[i];
  }

  inline Index CartesianProduct::startT(size_t i) const
  {
    assert(i < numberOfSubmanifolds() && "invalid index");
    return startIndexT_[i];
  }
}

#endif //_PGS_CARTESIAN_PRODUCT_H_

