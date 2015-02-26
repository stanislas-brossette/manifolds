#include <stdexcept>
#include <manifolds/Manifold.h>
#include <manifolds/pgs_assert.h>

namespace pgs
{
  Manifold::Manifold(Index dimension, Index tangentDimension, Index representationDimension)
    : dimension_(dimension)
    , tangentDim_(tangentDimension)
    , representationDim_(representationDimension)
    , lock_(false)
  {
    pgs_assert(dimension>=0 && "Negative dimension not accepted");
    pgs_assert(tangentDimension >= 0 && "Negative dimension not accepted");
    pgs_assert(representationDimension >= 0 && "Negative dimension not accepted");
    pgs_assert(dimension <= tangentDimension);
    pgs_assert(tangentDimension <= representationDimension);
  }
  
  Manifold::~Manifold()
  {
  }

  Point Manifold::createPoint() const
  {
    pgs_assert(isValid() || seeMessageAbove());
    lock();
    return Point(*this);
  }

  Point Manifold::createPoint(const ConstRefVec& val) const
  {
    pgs_assert(isValid() || seeMessageAbove());
    if (isValidInit(val))
    {
      lock();
      return Point(*this, val);
    }
    else
    {
      throw std::runtime_error("Bad Point Initialization");
    }
  }

  bool Manifold::isValidInit(const Eigen::VectorXd& val) const
  {
    pgs_assert(isValid() || seeMessageAbove());
    pgs_assert(val.size() == representationDim());
    return isValidInit_(val);
  }

  Point Manifold::getIdentity() const
  {
    pgs_assert(isValid() || seeMessageAbove());
    lock();
    Eigen::VectorXd id(representationDim_);
    setIdentity(id);
    return Point(*this, id);
  }

  Index Manifold::dim() const
  {
    pgs_assert(isValid() || seeMessageAbove());
    return dimension_;
  }

  Index Manifold::representationDim() const
  {
    pgs_assert(isValid() || seeMessageAbove());
    return representationDim_;
  }

  void Manifold::plus(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const
  {
    pgs_assert(isValid() || seeMessageAbove());
    pgs_assert(out.size() == representationDim_);
    pgs_assert(x.size() == representationDim_);
    pgs_assert(v.size() == dimension_);
    plus_(out, x, v);
  }

  void Manifold::minus(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const
  {
    pgs_assert(isValid() || seeMessageAbove());
    pgs_assert(out.size() == dimension_);
    pgs_assert(x.size() == representationDim_);
    pgs_assert(y.size() == representationDim_);
    minus_(out, x, y);
  }

  void Manifold::invMap(RefVec out, const ConstRefVec& x) const
  {
    pgs_assert(out.size() == dimension_);
    pgs_assert(x.size() == representationDim_);
    invMap_(out, x);
  }

  Eigen::MatrixXd Manifold::diffMap(const ConstRefVec& x) const
  {
    pgs_assert(isValid() || seeMessageAbove());
    pgs_assert(x.size() == representationDim_);
    return diffMap_(x);
  }

  void Manifold::applyDiffMap(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    pgs_assert(isValid() || seeMessageAbove());
    pgs_assert(in.cols() == representationDim_);
    pgs_assert(out.cols() == dimension_);
    pgs_assert(in.rows() == out.rows());
    pgs_assert(x.size() == representationDim_);
    applyDiffMap_(out, in, x);
  }

  Eigen::MatrixXd Manifold::diffInvMap(const ConstRefVec& x) const
  {
    pgs_assert(isValid() || seeMessageAbove());
    pgs_assert(x.size() == representationDim_);
    return diffInvMap_(x);
  }

  void Manifold::applyDiffInvMap(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x) const
  {
    pgs_assert(isValid() || seeMessageAbove());
    pgs_assert(out.cols() == representationDim_);
    pgs_assert(in.cols() == dimension_);
    pgs_assert(in.rows() == out.rows());
    pgs_assert(x.size() == representationDim_);
    applyDiffInvMap_(out, in, x);
  }

  void Manifold::applyTransport(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    pgs_assert(isValid() || seeMessageAbove());
    pgs_assert(in.rows() == dimension_);
    pgs_assert(out.rows() == dimension_);
    pgs_assert(in.cols() == out.cols());
    pgs_assert(x.size() == representationDim());
    pgs_assert(v.size() == dim());
    applyTransport_(out, in, x, v);
  }

  void Manifold::applyInvTransport(
      RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const
  {
    pgs_assert(isValid() || seeMessageAbove());
    pgs_assert(in.cols() == dimension_);
    pgs_assert(out.cols() == dimension_);
    pgs_assert(in.rows() == out.rows());
    pgs_assert(x.size() == representationDim());
    pgs_assert(v.size() == dim());
    applyInvTransport_(out, in, x, v);
  }

  void Manifold::setDimension(Index d)
  {
    pgs_assert(isValid() || seeMessageAbove());
    pgs_assert(d>=0 && "Negative dimension not accepted");
    testLock();
    dimension_ = d;
  }

  void Manifold::setRepresentationDimension(Index rd)
  {
    pgs_assert(isValid() || seeMessageAbove());
    pgs_assert(rd>=0 && "Negative dimension not accepted");
    testLock();
    representationDim_ = rd;
  }


  void Manifold::setIdentity(RefVec out) const
  {
    pgs_assert(isValid() || seeMessageAbove());
    pgs_assert(out.size() == static_cast<int> (representationDim_));
    setIdentity_(out);
  }


  void Manifold::lock() const
  {
    lock_ = true;
  }

  void Manifold::testLock() const
  {
    if (lock_)
      throw std::runtime_error("Either a point or a compound manifold is relying on this manifold, you can't modify it anymore.");
  }
}
