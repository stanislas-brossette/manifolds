#include <iostream>
#include <Eigen/Core>
#include <pgsolver/manifolds/defs.h>
#include <pgsolver/solver/Problem.h>

namespace pgs
{
  Problem::Problem(Manifold& manifold)
    :M_(manifold),
     x_(manifold.getIdentity()),
     z_(Eigen::VectorXd::Zero(manifold.dim())),
     phi_x_z_(x_+z_)
  {
  }

  Problem::Problem(Manifold& manifold, const Point& x)
    :M_(manifold),
     x_(x),
     z_(Eigen::VectorXd::Zero(manifold.dim())),
     phi_x_z_(x_+z_)
  {
    assert(&M_ == &(x.getManifold()) && "new point must belong to the problem's manifold");
  }

  void Problem::setX(const Point& x)
  {
    assert(&M_ == &(x.getManifold()) && "new point must belong to the problem's manifold");
    x_ = x;
    phi_x_z_ = x_ + z_;
    broadcastXIsNew();
  }

  void Problem::setZ(const ConstRefVec& z)
  {
    assert(z.size() == M_.dim() && "Wrong point increment value for this problem");
    z_ = z;
    phi_x_z_ = x_ + z_;
    broadcastZIsNew();
  }

  const Point& Problem::x() const
  {
    return x_;
  }

  const Point& Problem::phi_x_z() const
  {
    return phi_x_z_;
  }

  const Eigen::VectorXd& Problem::z() const
  {
    return z_;
  }

  const Manifold& Problem::M() const
  {
    return M_;
  }

  void Problem::broadcastXIsNew()
  {
    //Do Nothing
  }

  void Problem::broadcastZIsNew()
  {
    //Do Nothing
  }

  Index Problem::linCstrDim() const
  {
    Index tot = 0;
    for(size_t i = 0; i < numberOfCstr(); ++i)
      tot += nonLinCstrDim(i);
    return tot;
  }

  Index Problem::nonLinCstrDim() const
  {
    Index tot = 0;
    for(size_t i = 0; i < numberOfCstr(); ++i)
      tot += nonLinCstrDim(i);
    return tot;
  }

  void Problem::evalLinCstr(RefVec out) const
  {
    Index totLinDim = linCstrDim();
    assert(out.size() == totLinDim && "wrong total size for linear cstr");
    Index startIndex = 0;
    for(size_t i = 0; i < numberOfCstr(); ++i)
    {
      evalLinCstr(out.segment(startIndex,linCstrDim(i)), i);
      startIndex += linCstrDim(i);
    }
  }
  void Problem::evalLinCstrGrad(RefMat out) const
  {
    assert(out.rows() == linCstrDim() && "wrong total size for linear cstr");
    assert(out.cols() == M().dim() && "Wrong cols size");
    Index startIndex = 0;
    for(size_t i = 0; i < numberOfCstr(); ++i)
    {
      evalLinCstrGrad(out.middleRows(startIndex,linCstrDim(i)), i);
      startIndex += linCstrDim(i);
    }
  }
  void Problem::getLinCstrLB(RefVec out) const
  {
    Index totLinDim = linCstrDim();
    assert(out.size() == totLinDim && "wrong total size for linear cstr");
    Index startIndex = 0;
    for(size_t i = 0; i < numberOfCstr(); ++i)
    {
      getLinCstrLB(out.segment(startIndex,linCstrDim(i)), i);
      startIndex += linCstrDim(i);
    }
  }
  void Problem::getLinCstrUB(RefVec out) const
  {
    Index totLinDim = linCstrDim();
    assert(out.size() == totLinDim && "wrong total size for linear cstr");
    Index startIndex = 0;
    for(size_t i = 0; i < numberOfCstr(); ++i)
    {
      getLinCstrUB(out.segment(startIndex,linCstrDim(i)), i);
      startIndex += linCstrDim(i);
    }
  }

  void Problem::evalNonLinCstr(RefVec out) const
  {
    Index totNonLinDim = nonLinCstrDim();
    assert(out.size() == totNonLinDim && "wrong total size for Nonlinear cstr");
    Index startIndex = 0;
    for(size_t i = 0; i < numberOfCstr(); ++i)
    {
      evalNonLinCstr(out.segment(startIndex,nonLinCstrDim(i)), i);
      startIndex += nonLinCstrDim(i);
    }
  }
  void Problem::evalNonLinCstrGrad(RefMat out) const
  {
    assert(out.rows() == nonLinCstrDim() && "wrong total size for Nonlinear cstr");
    assert(out.cols() == M().dim() && "Wrong cols size");
    Index startIndex = 0;
    for(size_t i = 0; i < numberOfCstr(); ++i)
    {
      evalNonLinCstrGrad(out.middleRows(startIndex,nonLinCstrDim(i)), i);
      startIndex += nonLinCstrDim(i);
    }
  }
  void Problem::getNonLinCstrLB(RefVec out) const
  {
    Index totNonLinDim = nonLinCstrDim();
    assert(out.size() == totNonLinDim && "wrong total size for nonlinear cstr");
    Index startIndex = 0;
    for(size_t i = 0; i < numberOfCstr(); ++i)
    {
      getNonLinCstrLB(out.segment(startIndex,nonLinCstrDim(i)), i);
      startIndex += nonLinCstrDim(i);
    }
  }
  void Problem::getNonLinCstrUB(RefVec out) const
  {
    Index totNonLinDim = nonLinCstrDim();
    assert(out.size() == totNonLinDim && "wrong total size for nonlinear cstr");
    Index startIndex = 0;
    for(size_t i = 0; i < numberOfCstr(); ++i)
    {
      getNonLinCstrUB(out.segment(startIndex,nonLinCstrDim(i)), i);
      startIndex += nonLinCstrDim(i);
    }
  }
  
  void Problem::printState() const
  {
    double obj = 0.0;
    Eigen::MatrixXd objGrad(1,M().dim());
    Eigen::VectorXd tangentLB(M().dim());
    Eigen::VectorXd tangentUB(M().dim());

    Index dimLin = linCstrDim();
    Eigen::VectorXd linCstr(dimLin);
    Eigen::MatrixXd linCstrDiff(dimLin,M().dim());
    Eigen::VectorXd linCstrLB(dimLin);
    Eigen::VectorXd linCstrUB(dimLin);

    Index dimNonLin = nonLinCstrDim();
    Eigen::VectorXd NonLinCstr(dimNonLin);
    Eigen::MatrixXd NonLinCstrDiff(dimNonLin,M().dim());
    Eigen::VectorXd NonLinCstrLB(dimNonLin);
    Eigen::VectorXd NonLinCstrUB(dimNonLin);

    std::cout << "Current status of the problem" << std::endl;
    std::cout << "x=" << x() << std::endl;
    std::cout << "z=" << z().transpose() << std::endl;
    std::cout << "x + z = " << x()+z() << std::endl;

    std::cout << std::endl << "Objective Function:" << std::endl;
    evalObj(obj);
    evalObjGrad(objGrad);
    std::cout << "f(phi_x(z))=" << obj << std::endl;
    std::cout << "grad_z(f(phi_x(z))=" << objGrad << std::endl;

    std::cout << std::endl << "Bounds:" << std::endl;
    getTangentLB(tangentLB);
    getTangentUB(tangentUB);
    std::cout << "[" << tangentLB.transpose() << "] < z < [" << tangentUB.transpose() << "]" << std::endl;

    std::cout << std::endl << "Linear Constraints:" << std::endl;
    evalLinCstr(linCstr);
    evalLinCstrGrad(linCstrDiff);
    getLinCstrLB(linCstrLB);
    getLinCstrUB(linCstrUB);
    std::cout << linCstrLB << " < linCstr = " << linCstr << " < " << linCstrUB << std::endl;
    std::cout << "gradlinCstr = " << linCstrDiff << std::endl;

    std::cout << std::endl << "NonLinear Constraints:" << std::endl;
    evalNonLinCstr(NonLinCstr);
    evalNonLinCstrGrad(NonLinCstrDiff);
    getNonLinCstrLB(NonLinCstrLB);
    getNonLinCstrUB(NonLinCstrUB);
    std::cout << NonLinCstrLB << " < NonLinCstr = " << NonLinCstr << " < " << NonLinCstrUB << std::endl;
    std::cout << "gradNonLinCstr = " << NonLinCstrDiff << std::endl;
  }
}
