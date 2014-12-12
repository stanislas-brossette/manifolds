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
      /// \brief Default Constructor
      Problem(const Manifold& manifold);
      
      /// \brief Get Upper Bounds on the variable
      virtual void getUB(RefVec out) = 0;
      /// \brief Get Lower Bounds on the variable
      virtual void getLB(RefVec out) = 0;

      /// \brief Get Constraints Lower Bounds
      virtual void getCstrLB(RefVec out) = 0;
      /// \brief Get Constraints Upper Bounds
      virtual void getCstrUB(RefVec out) = 0;

      /// \brief Evaluate Objective Function at point
      //\f$\phi_x^{\mathcal{M}}(p)\f$
      virtual void evalObj(RefVec out) = 0;
      /// \brief Evaluate Gradient of Objective Function at point
      //\f$x\f$
      virtual void evalObjGrad(RefVec out) = 0;

      /// \brief Evaluate Linear Constraints Index i at point
      //\f$\phi_x^{\mathcal{M}}(p)\f$
      virtual void evalLinCstr(RefVec out, Index i) = 0;
      /// \brief Evaluate Gradient of Linear Constraints Index i\n
      /// They are constants
      virtual void evalLinCstrGrad(RefVec out, Index i) = 0;

      /// \brief Evaluate NonLinear Constraints Index i at point
      //\f$\phi_x^{\mathcal{M}}(p)\f$
      virtual void evalNonLinCstr(RefVec out, Index i) = 0;
      /// \brief Evaluate Gradient of NonLinear Constraints Index i at point
      //\f$x\f$
      virtual void evalNonLinCstrGrad(RefVec out, Index i) = 0;

    private:
      /// \brief Manifold on which the problem is defined
      const Manifold& M_;
      /// \brief Current zero point of the map
      Point x_;
      /// \brief Current increment
      Eigen::VectorXd p_;
  };
}

#endif //_PGS_PROBLEM_H_
