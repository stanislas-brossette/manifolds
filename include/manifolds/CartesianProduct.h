#ifndef _MANIFOLDS_CARTESIAN_PRODUCT_H_
#define _MANIFOLDS_CARTESIAN_PRODUCT_H_

#include <vector>
#include <manifolds/defs.h>
#include <manifolds/Manifold.h>

namespace pgs
{
  /// \brief Manifold representing the cartesian product of several submanifolds
  class PGS_API CartesianProduct : public Manifold
  {
  public:
    /// \brief Default constructor
    CartesianProduct();

    /// \brief Constructor of the manifold composed of \f$ m1\times m2\f$
    CartesianProduct(const Manifold& m1, const Manifold& m2);

    /// \brief Adds manifold m to the current composed manifold\n
    /// This method cannot be executed if the manifold is locked
    CartesianProduct& multiply(const Manifold& m);

    virtual size_t numberOfSubmanifolds() const;
    virtual const Manifold& operator()(size_t i) const;

    virtual std::string toString(const ConstRefVec& val, const std::string& prefix = "") const;

  protected:
    virtual bool isValidInit_(const Eigen::VectorXd& val) const;

    virtual Index startR(size_t i) const;
    virtual Index startT(size_t i) const;

    virtual void plus_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void minus_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const;
    virtual void invMap_(RefVec out, const ConstRefVec& x) const;
    virtual void setIdentity_(RefVec out) const;
    virtual Eigen::MatrixXd diffMap_(const ConstRefVec& x) const;
    virtual void applyDiffMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    virtual Eigen::MatrixXd diffInvMap_(const ConstRefVec& x) const;
    virtual void applyDiffInvMap_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    virtual void applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const;

  private:
    /// \brief List of pointers on all the manifolds in the cartesian product
    std::vector<const Manifold* > submanifolds_;

    /// \brief List of start index of submanifolds in a vector of the
    /// tangent space
    std::vector<Index> startIndexT_;

    /// \brief List of start index of submanifolds in a vector of the
    /// representation space
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

#endif //_MANIFOLDS_CARTESIAN_PRODUCT_H_

