#ifndef _PGS_PROBLEM_H_
#define _PGS_PROBLEM_H_

#include <iostream>
#include <Eigen/Core>

#include <pgsolver/manifolds/defs.h>
#include <pgsolver/manifolds/Point.h>
#include <pgsolver/manifolds/Manifold.h>

namespace pgs
{
  class Problem
  {
    public:
      /// \brief Constructor taking a manifold a setting x and z to zero
      Problem(Manifold& manifold);
      /// \brief Constructor that sets x
      Problem(Manifold& manifold, const Point& x);
      
      /// \brief Sets new value for x_
      void setX(const Point& x);
      /// \brief Sets new value for z_
      void setZ(const ConstRefVec& z);

      /// \brief Gets current x_ value
      const Point& x() const;
      /// \brief Gets current z_ value
      const Eigen::VectorXd& z() const;
      /// \brief Gets current manifold M_ value
      const Manifold& M() const;
      /// \brief Gets current value of \f$\phi_x(z)
      const Point& phi_x_z() const;
      
      /// \brief Get Lower Bounds on the variable\n
      /// Correspond to the bounds on z for a given x.
      virtual void getTangentLB(RefVec out) const = 0;
      /// \brief Get Upper Bounds on the variable
      /// Correspond to the bounds on z for a given x.
      virtual void getTangentUB(RefVec out) const = 0;

      /// \brief Evaluate Objective Function at point
      //\f$\phi_x^{\mathcal{M}}(z)\f$
      virtual void evalObj(double& out) const = 0;
      /// \brief Evaluate Gradient of Objective Function at point
      //\f$x\f$
      virtual void evalObjGrad(RefMat out) const = 0;

      /// \brief Evaluate Linear Constraints Index i at point
      /// \f$\phi_x^{\mathcal{M}}(z)\f$
      virtual void evalLinCstr(RefVec out, size_t i) const = 0;
      /// \brief Evaluate Gradient of Linear Constraints Index i\n
      /// They are constants
      virtual void evalLinCstrGrad(RefMat out, size_t i) const = 0;
      /// \brief Get Constraints Lower Bounds
      virtual void getLinCstrLB(RefVec out, size_t i) const = 0;
      /// \brief Get Constraints Upper Bounds
      virtual void getLinCstrUB(RefMat out, size_t i) const = 0;

      /// \brief Evaluate NonLinear Constraints Index i at point
      //\f$\phi_x^{\mathcal{M}}(z)\f$
      virtual void evalNonLinCstr(RefVec out, size_t i) const = 0;
      /// \brief Evaluate Gradient of NonLinear Constraints Index i at point x
      virtual void evalNonLinCstrGrad(RefMat out, size_t i) const = 0;
      /// \brief Get Constraints Lower Bounds
      virtual void getNonLinCstrLB(RefVec out, size_t i) const = 0;
      /// \brief Get Constraints Upper Bounds
      virtual void getNonLinCstrUB(RefVec out, size_t i) const = 0;

    protected:
      /// \brief Updates the problem for a new value of X
      virtual void broadcastXIsNew();
      /// \brief Updates the problem for a new value of Z
      virtual void broadcastZIsNew();

    protected:
      /// \brief Manifold on which the problem is defined
      Manifold& M_;

    private:
      /// \brief Current zero point of the map
      Point x_;
      /// \brief Current tangent vector \f$ z \in T_x\mathcal{M}\f$
      Eigen::VectorXd z_;
      /// \brief Current Point x+z
      Point phi_x_z_;
  };
}

#endif //_PGS_PROBLEM_H_
