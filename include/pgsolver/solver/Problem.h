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
      /// \brief Evaluate Jacobian of Objective Function at point
      //\f$x\f$
      virtual void evalObjDiff(RefMat out) const = 0;

      /// \brief Evaluate Linear Constraints Index i at point
      /// \f$\phi_x^{\mathcal{M}}(z)\f$
      virtual void evalLinCstr(RefMat out, size_t i) const = 0;
      /// \brief Evaluate Jacobian of Linear Constraints Index i\n
      /// They are constants
      virtual void evalLinCstrDiff(RefMat out, size_t i) const = 0;
      /// \brief Get Constraints Lower Bounds
      virtual void getLinCstrLB(RefMat out, size_t i) const = 0;
      /// \brief Get Constraints Upper Bounds
      virtual void getLinCstrUB(RefMat out, size_t i) const = 0;

      /// \brief Evaluate NonLinear Constraints Index i at point
      //\f$\phi_x^{\mathcal{M}}(z)\f$
      virtual void evalNonLinCstr(RefMat out, size_t i) const = 0;
      /// \brief Evaluate Jacobian of NonLinear Constraints Index i at point x
      virtual void evalNonLinCstrDiff(RefMat out, size_t i) const = 0;
      /// \brief Get Constraints Lower Bounds
      virtual void getNonLinCstrLB(RefMat out, size_t i) const = 0;
      /// \brief Get Constraints Upper Bounds
      virtual void getNonLinCstrUB(RefMat out, size_t i) const = 0;

      /// \brief Evaluate All Linear Constraints at point
      /// \f$\phi_x^{\mathcal{M}}(z)\f$
      virtual void evalLinCstr(RefVec out) const;
      /// \brief Evaluate All Jacobian of Linear Constraints\n
      /// They are constants
      virtual void evalLinCstrDiff(RefMat out) const;
      /// \brief Get All Constraints Lower Bounds
      virtual void getLinCstrLB(RefVec out) const;
      /// \brief Get All Constraints Upper Bounds
      virtual void getLinCstrUB(RefVec out) const;

      /// \brief Evaluate All NonLinear Constraints at point
      //\f$\phi_x^{\mathcal{M}}(z)\f$
      virtual void evalNonLinCstr(RefVec out) const;
      /// \brief Evaluate Jacobian of All NonLinear Constraints at point x
      virtual void evalNonLinCstrDiff(RefMat out) const;
      /// \brief Get All Constraints Lower Bounds
      virtual void getNonLinCstrLB(RefVec out) const;
      /// \brief Get All Constraints Upper Bounds
      virtual void getNonLinCstrUB(RefVec out) const;

      /// \brief Function displaying the current state of the problem
      virtual void printState() const;
      /// \brief Total Number of Consraints
      virtual size_t numberOfCstr() const = 0;
      /// \brief Dimension of linear part of Consraints i
      virtual Index linCstrDim(size_t i) const = 0;
      /// \brief Computes the sum of Dimension of linear part of all Consraints
      virtual Index linCstrDim() const;
      /// \brief Dimension of nonlinear part of Consraints i
      virtual Index nonLinCstrDim(size_t i) const = 0;
      /// \brief Computes the sum of Dimension of nonlinear part of all Consraints
      virtual Index nonLinCstrDim() const;

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
      /// \brief Current Point \f$ x \oplus z = \phi_x(z) \f$ 
      Point phi_x_z_;
  };
}

#endif //_PGS_PROBLEM_H_
