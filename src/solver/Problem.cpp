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
    std::cout << "Problem constructor"<< std::endl;
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

}
