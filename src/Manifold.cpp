#include <stdexcept>
#include <pgsolver/Manifold.h>

namespace pgs
{
  Manifold::Manifold(Index dimension, Index representationDimension)
    : dimension_(dimension)
    , representationDim_(representationDimension)
    , lock_(false)
  {
    assert(dimension>=0 && "Negative dimension not accepted");
    assert(representationDimension>=0 && "Negative dimension not accepted");
  }

  Point Manifold::createPoint() const
  {
    lock();
    incrementRefCounter();
    return Point(*this);
  }

  Point Manifold::createPoint(const Eigen::VectorXd& val) const
  {
    if(isValidInit(val))
    {
      lock();
      incrementRefCounter();
      return Point(*this, val);
    }
    else
    {
      throw std::runtime_error("Bad Point Initialization"); 
    }
  }

  Point Manifold::getIdentity() const
  {
    lock();
    incrementRefCounter();
    Eigen::VectorXd id(representationDim_);
    setIdentity(id);
    return Point(*this, id);
  }

  Index Manifold::dim() const
  {
    return dimension_;
  }

  Index Manifold::representationDim() const
  {
    return representationDim_;
  }

  void Manifold::plus(RefVec out, ConstRefVec& x, ConstRefVec& v) const
  {
    assert(out.size() ==  representationDim_);
    assert(x.size() == representationDim_);
    assert(v.size() == dimension_);
    plus_(out, x, v);
  }

  void Manifold::minus(RefVec out, ConstRefVec& x, ConstRefVec& y) const
  {
    assert(out.size() == dimension_);
    assert(x.size() == representationDim_);
    assert(y.size() == representationDim_);
    minus_(out, x, y);
  }

  Eigen::MatrixXd Manifold::diffMap(ConstRefVec& x) const
  {
    assert(x.size() == representationDim_);
    return diffMap_(x);
  }
  
  void Manifold::applyDiffMap(RefMat inOut, ConstRefVec& x) const
  {
    assert(inOut.rows() == representationDim_);
    assert(inOut.cols() == dimension_);
    assert(x.size() == representationDim_);
    applyDiffMap_(inOut, x);
  }

  void Manifold::setDimension(Index d)
  {
    assert(d>0 && "Negative dimension not accepted");
    testLock();
    dimension_ = d;
  }

  void Manifold::setRepresentationDimension(Index rd)
  {
    assert(rd>0 && "Negative dimension not accepted");
    testLock();
    representationDim_ = rd;
  }


  void Manifold::setIdentity(RefVec out) const
  {
    assert(out.size() == static_cast<int> (representationDim_));
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
