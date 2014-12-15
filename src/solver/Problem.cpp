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
  
  void Problem::printState() const
  {
    double obj = 0.0;
    Eigen::MatrixXd objGrad(1,M().dim());
    Eigen::VectorXd tangentLB(M().dim());
    Eigen::VectorXd tangentUB(M().dim());

    Eigen::VectorXd linCstr(nbLin);
    Eigen::MatrixXd linCstrDiff(nbLin,M().dim());
    Eigen::VectorXd linCstrLB(nbLin);
    Eigen::VectorXd linCstrUB(nbLin);

    Eigen::VectorXd NonLinCstr(nbNonLin);
    Eigen::MatrixXd NonLinCstrDiff(nbNonLin,M().dim());
    Eigen::VectorXd NonLinCstrLB(nbNonLin);
    Eigen::VectorXd NonLinCstrUB(nbNonLin);

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
    evalLinCstr(linCstr, 0);
    evalLinCstrGrad(linCstrDiff, 0);
    getLinCstrLB(linCstrLB, 0);
    getLinCstrUB(linCstrUB, 0);
    std::cout << linCstrLB << " < linCstr = " << linCstr << " < " << linCstrUB << std::endl;
    std::cout << "gradlinCstr = " << linCstrDiff << std::endl;

    std::cout << std::endl << "NonLinear Constraints:" << std::endl;
    evalNonLinCstr(NonLinCstr, 0);
    evalNonLinCstrGrad(NonLinCstrDiff, 0);
    getNonLinCstrLB(NonLinCstrLB, 0);
    getNonLinCstrUB(NonLinCstrUB, 0);
    std::cout << NonLinCstrLB << " < NonLinCstr = " << NonLinCstr << " < " << NonLinCstrUB << std::endl;
    std::cout << "gradNonLinCstr = " << NonLinCstrDiff << std::endl;
  }

}
