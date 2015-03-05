#ifndef _MANIFOLDS_CARTESIAN_PRODUCT_H_
#define _MANIFOLDS_CARTESIAN_PRODUCT_H_

#include <vector>
#include <manifolds/defs.h>
#include <manifolds/Manifold.h>

namespace pgs
{
  /// \brief Manifold representing the cartesian product of several submanifolds
  class MANIFOLDS_API CartesianProduct : public Manifold
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
    virtual bool isInM_(const Eigen::VectorXd& val) const;

    virtual Index startR(size_t i) const;
    virtual Index startT(size_t i) const;

    virtual void retractation_(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void pseudoLog_(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const;
    virtual void invMap_(RefVec out, const ConstRefVec& x) const;
    virtual void setZero_(RefVec out) const;
    virtual Eigen::MatrixXd diffRetractation_(const ConstRefVec& x) const;
    virtual void applyDiffRetractation_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    virtual Eigen::MatrixXd diffPseudoLog0_(const ConstRefVec& x) const;
    virtual void applyDiffPseudoLog0_(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;
    virtual void applyTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void applyInvTransport_(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const;
    
    virtual void tangentConstraint_(RefMat out, const ConstRefVec& x) const;
    virtual bool isInTxM_(const ConstRefVec& x, const ConstRefVec& v) const;
    virtual void forceOnTxM_(RefVec out, const ConstRefVec& in, const ConstRefVec& x) const;

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
    pgs_assert(i < numberOfSubmanifolds() && "invalid index");
    return startIndexR_[i];
  }

  inline Index CartesianProduct::startT(size_t i) const
  {
    pgs_assert(i < numberOfSubmanifolds() && "invalid index");
    return startIndexT_[i];
  }
}

#endif //_MANIFOLDS_CARTESIAN_PRODUCT_H_

