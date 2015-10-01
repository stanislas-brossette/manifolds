// Copyright (c) 2015 CNRS
// Authors: Stanislas Brossette, Adrien Escande

// This file is part of manifolds
// manifolds is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.

// manifolds is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// manifolds. If not, see
// <http://www.gnu.org/licenses/>.

#pragma once

#include <iostream>

#include <Eigen/Core>

#include <manifolds/defs.h>
#include <manifolds/view.h>
//#include <manifolds/RefCounter.h>
//#include <manifolds/Point.h>
//#include <manifolds/ValidManifold.h>
#include <manifolds/Manifold_Base.h>

namespace mnf
{
  class Point;
  /// \brief The Manifold Class represents a manifold. It contains the implementations of
  /// the basic operations on it, like external addition, internal substraction,
  /// Translation from tangent space to representation space and back, derivatives
  /// etc...\n
  /// Manifold derives from RefCounter to keep track of how many points depend on it\n
  /// A manifold is characterized by its tangent space and its representation
  /// space and by the map to go from a vector of the tangent space to a point
  /// of the manifold.
  /// We denote the a manifold as \f$ \mathcal{M}\f$, its tangent space at point
  /// \f$ x \f$ is \f$ T_x^\mathcal{M} \f$ seen as subspace of \f$\mathbb{R}^t \f$ and its
  /// representation space is denoted \f$\mathbb{E}\f$\n
  /// The map function is \f$ \phi:\mathbb{M},T^\mathbb{M}\to\mathbb{M}
  /// \f$ and the map function on a point \f$x\f$ is \f$ \phi_x:T_x^\mathbb{M}\to\mathbb{M}\f$
  class MANIFOLDS_API ConstManifold
  {
  protected:
    ConstManifold(std::shared_ptr<Manifold_Base> m);
  public:
    ConstManifold(std::shared_ptr<const Manifold_Base> m);

    std::shared_ptr<const Manifold_Base> ptr() const;

  public:
    /// \brief The destructor
    virtual ~ConstManifold();

    /// \brief CreatePoint allows to create a point that belongs to a manifold and that
    /// behaves according to the manifolds operations
    Point createPoint() const;

    /// \brief Overload of createPoint to set the value of the point
    Point createPoint(const ConstRefVec& val) const;

    /// \brief Creates a point that represents the identity wrt the addition operation
    /// defined in this manifold (aka the zero)
    Point getZero() const;

    /// \brief Creates a random point on this manifold (aka the retractation from point 0
    /// of a random vector of the tangent space)
    Point createRandomPoint(double coeff = 1.0) const;
    void createRandomPoint(RefVec out, double coeff = 1.0) const;

    /// \brief Checks that the value val described in the representation space
    /// is an element of the manifold
    virtual bool isInM(const Eigen::VectorXd& val, double prec = 1e-8) const;

    /// \brief finds the closest point to \a in on \f$ \mathbb{M} \f$.
    virtual void forceOnM(RefVec out, const ConstRefVec& in) const;

    /// \brief Returns the dimension of the Manifold
    Index dim() const;

    /// \brief Returns the number \a t of parameters used to represent a tangent vector
    /// i.e. T_x\mathbb{M} is seen as a subspace of \f$ \mathbb{R}^t \f$
    Index tangentDim() const;

    /// \brief Returns the dimension of the representation space of the manifold
    Index representationDim() const;

    /// \brief Returns True if the manifold is an elementary manyfold, false otherwise
    virtual bool isElementary() const;

    /// \brief Displays a description of the manifold
    virtual void display(const std::string& prefix = "") const;

    /// \brief Returns the number of submanifolds that compose the manifold
    virtual size_t numberOfSubmanifolds() const;

    /// \brief Returns a pointer on the submanifold of index i if it exists
    /// Only useful with composed Manifolds
    virtual ConstManifold operator()(size_t i) const;

    //view
    /// \brief Returns a view of vector val as seen as an element of submanifold i.
    /// Template D indicates if we considere the tangent space or representation
    /// space of submanifold i
    template<int D> Segment getView(RefVec val, size_t i) const;

    /// \brief Const version of getView
    template<int D> ConstSegment getConstView(const ConstRefVec& val, size_t i) const;

    template<int Dr, int Dc> typename ViewReturnType<Dr, Dc>::Type getView(RefMat val, size_t i) const;
    template<int Dr, int Dc> typename ConstViewReturnType<Dr, Dc>::Type getConstView(const ConstRefMat& val, size_t i) const;
    template<int Dr, int Dc> typename ViewReturnType<Dr, Dc>::Type getView(RefMat val, size_t ir, size_t ic) const;
    template<int Dr, int Dc> typename ConstViewReturnType<Dr, Dc>::Type getConstView(const ConstRefMat& val, size_t ir, size_t ic) const;

    /// \brief returns a boolean indicating whether this manifold is currently locked
    bool isLocked() const;

    /// \brief Converts val to string for pretty printing
    virtual std::string toString(const ConstRefVec& val, const std::string& prefix = "", const Eigen::IOFormat& fmt = mnf::defaultFormat) const;

    //map operations
    //The following operations are public and make call to their private version
    //that actually contain the implementation

    /// \brief Sets vector out to the zero of this manifold
    void setZero(RefVec out) const;

    /// \brief External addition operation \f$ out = x \oplus v = \phi_x(v) \f$
    /// \param out output reference on element of the manifold
    /// \f$out\in\mathbb{M}\f$
    /// \param x element of the manifold \f$x\in\mathbb{M}\f$
    /// \param v element of the tangent space of the manifold
    /// \f$v\in T_x^\mathcal{M}\f$
    void retractation(RefVec out, const ConstRefVec& x, const ConstRefVec& v) const;

    /// \brief PseudoLog operation
    /// \f$ out = {Log}_x(y) \f$
    /// \param out output reference on element of the tangent space of the
    /// manifold \f$out\in\mathbb{M}\f$
    /// \param x element of the manifold \f$x\in\mathbb{M}\f$
    /// \param y element of the manifold \f$y\in\mathbb{M}\f$
    void pseudoLog(RefVec out, const ConstRefVec& x, const ConstRefVec& y) const;

    /// \brief computes the inverse of a point of the manifold through its map
    /// \f$ out = x \ominus 0 = \phi_0^{-1}(x) \f$
    /// \param out output reference on element of the tangent space of the
    /// manifold\f$out\in\mathbb{M}\f$
    /// \param x element of the manifold\f$x\in\mathbb{M}\f$
    void pseudoLog0(RefVec out, const ConstRefVec& x) const;

    /// \brief Computes the Jacobian matrix of the map function
    /// \f$\frac{\partial\phi_x}{\partial v}(0)\f$
    /// \param x element of manifold \f$x\in\mathbb{M}\f$
    /// \return Matrix representing \f$\frac{\partial\phi_x}{\partial v}(0)\f$
    Eigen::MatrixXd diffRetractation(const ConstRefVec& x) const;

    /// \brief Computes the product of a matrix in with the jacobian matrix of
    /// the map on point x.\n \f$ out = in*\frac{\partial\phi_x}{\partial v}(0)\f$
    /// \param out result of the operation
    /// \param in matrix to which the operation is applied,
    /// each row of in must be a vector of the ambient space
    /// \param x point of the manifold on which the map is taken
    void applyDiffRetractation(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;

    /// \brief Computes the Jacobian matrix of the pseudoLog0 function
    /// \f$\frac{\partial\phi^{-1}_0}{\partial x}(x)\f$
    /// \param x element of manifold \f$x\in\mathbb{M}\f$
    /// \return Matrix representing \f$\frac{\partial\phi^{-1}_0}{\partial x}(x)\f$
    Eigen::MatrixXd diffPseudoLog0(const ConstRefVec& x) const;

    /// \brief Computes the product of a matrix in with the jacobian matrix of
    /// the inverse map on point x.\n \f$ out = in*\frac{\partial\phi^{-1}_0}{\partial x}(x)\f$
    /// \param out result of the operation
    /// \param in matrix to which the operation is applied
    /// \param x point of the manifold on which the map is taken
    void applyDiffPseudoLog0(RefMat out, const ConstRefMat& in, const ConstRefVec& x) const;

    /// \brief applies a transport operation from point \f$x\in\mathcal{M}\f$ of
    /// direction \f$v\in T_x^\mathcal{M}\f$ on matrix in
    /// \param out result of the operation
    /// \param in matrix to which the operation is applied
    /// \param x element of manifold \f$x\in\mathbb{M}\f$
    /// \param v element of the tangent space of the manifold
    void applyTransport(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const;

    /// \brief applies an inverse transport operation from point \f$x\in\mathcal{M}\f$ of
    /// direction \f$v\in T_x^\mathcal{M}\f$ on matrix in
    /// \param out result of the operation
    /// \param in matrix to which the operation is applied
    /// \param x element of manifold \f$x\in\mathbb{M}\f$
    /// \param v element of the tangent space of the manifold
    void applyInvTransport(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const;

    /// \brief applies an inverse transport operation from point \f$x\in\mathcal{M}\f$ on the right side of a matrix
    /// direction \f$v\in T_x^\mathcal{M}\f$ on matrix in
    /// \param out result of the operation
    /// \param in matrix to which the operation is applied
    /// \param x element of manifold \f$x\in\mathbb{M}\f$
    /// \param v element of the tangent space of the manifold
    void applyInvTransportOnTheRight(RefMat out, const ConstRefMat& in, const ConstRefVec& x, const ConstRefVec& v) const;

    /// \brief returns the constraint matrix to test that a vector \f$v \in \mathbb{R}^t \f$
    /// belongs to \f$ T_x \mathbb{M} \f$. This is the case if \f$out v = 0\f$.
    /// Returns a (\a t - \a d) by \a t matrix.
    void tangentConstraint(RefMat out, const ConstRefVec& x) const;

    /// \brief checks if \a v belongs to \f$ T_x \mathbb{M} \f$.
    /// If \a t = \a d, return always true
    bool isInTxM(const ConstRefVec& x, const ConstRefVec& v, const double& prec = 1e-8) const;

    /// \brief finds the closest vector to \a in on \f$ T_x \mathbb{M} \f$.
    /// If \a t = \a d, \a out = \a in.
    void forceOnTxM(RefVec out, const ConstRefVec& in, const ConstRefVec&x) const;

    /// @brief Computes the Identity on TxM. Such that Id.v=v and every column is in TxM
    void getIdentityOnTxM(RefMat out, const ConstRefVec& x) const;

    /// \brief Locks the manifold\n Once a point of the manifold is created, the
    /// manifold cannot be modified anymore.
    void lock() const;

    /// \brief name_ getter
    const std::string& name() const;

    /// \brief fills the \a out vector with the limits of validity of the tangent map
    void limitMap(RefVec out) const;

    /// \brief Gets the typical Magnitude of the manifold on its tangent space.
    /// it is typically used for scaling of the variables
    void getTypicalMagnitude(RefVec out) const;
    Eigen::VectorXd getTypicalMagnitude() const;

    /// \brief Gets the trust Magnitude of the manifold on its tangent space.
    /// it is typically used for scaling of the trust regions
    void getTrustMagnitude(RefVec out) const;
    Eigen::VectorXd getTrustMagnitude() const;

    /// \brief returns a value unique to each instance of a manifold
    long getInstanceId() const;

    /// \brief Compares this manifold to another one using the identifier
    /// unique to each class implemented in getTypeId().
    virtual bool isSameType(const ConstManifold& other) const;

    /// \brief returns an id that should be unique to each manifold class.
    /// Use utils::hash::computeHash("some_string") to generate an ID at
    /// compile-time.
    virtual long getTypeId() const;

  protected:
    std::shared_ptr<Manifold_Base> ptr_;

  };
  class MANIFOLDS_API Manifold : public ConstManifold
  {
  public:
    Manifold(std::shared_ptr<Manifold_Base> m);
    std::shared_ptr<Manifold_Base> getNonConstPtr();
    /// \brief name_ setter
    std::string& name();
    void setTypicalMagnitude(const double& magnitude);
    void setTypicalMagnitude(const ConstRefVec& magnitude);
    void setTrustMagnitude(const double& magnitude);
    void setTrustMagnitude(const ConstRefVec& magnitude);
  };

  template<int D>
  inline Segment ConstManifold::getView(RefVec val, size_t i) const
  {
    return ptr_->getView<D>(val, i);
  }

  template<int D>
  inline ConstSegment ConstManifold::getConstView(const ConstRefVec& val, size_t i) const
  {
    return ptr_->getConstView<D>(val, i);
  }

  template<int Dr, int Dc>
  inline typename ViewReturnType<Dr, Dc>::Type ConstManifold::getView(RefMat val, size_t i) const
  {
    return ptr_->getView<Dr, Dc>(val, i);
  }

  template<int Dr, int Dc>
  inline typename ConstViewReturnType<Dr, Dc>::Type ConstManifold::getConstView(const ConstRefMat& val, size_t i) const
  {
    return ptr_->getConstView<Dr, Dc>(val, i);
  }

  template<int Dr, int Dc>
  inline typename ViewReturnType<Dr, Dc>::Type ConstManifold::getView(RefMat val, size_t ir, size_t ic) const
  {
    return ptr_->getView<Dr, Dc>(val, ir, ic);
  }

  template<int Dr, int Dc>
  inline typename ConstViewReturnType<Dr, Dc>::Type ConstManifold::getConstView(const ConstRefMat& val, size_t ir, size_t ic) const
  {
    return ptr_->getConstView<Dr, Dc>(val, ir, ic);
  }

  template<>
  inline typename ViewReturnType<F, R>::Type ConstManifold::getView<F, R>(RefMat val, size_t ir, size_t ic) const
  {
    return ptr_->getView<F, R>(val, ir, ic);
  }

  template<>
  inline typename ConstViewReturnType<F, R>::Type ConstManifold::getConstView<F, R>(const ConstRefMat& val, size_t ir, size_t ic) const
  {
    return ptr_->getConstView<F, R>(val, ir, ic);
  }

  template<>
  inline typename ViewReturnType<F, T>::Type ConstManifold::getView<F, T>(RefMat val, size_t ir, size_t ic) const
  {
    return ptr_->getView<F, T>(val, ir, ic);
  }

  template<>
  inline typename ConstViewReturnType<F, T>::Type ConstManifold::getConstView<F, T>(const ConstRefMat& val, size_t ir, size_t ic) const
  {
    return ptr_->getConstView<F, T>(val, ir, ic);
  }

  template<>
  inline typename ViewReturnType<R, F>::Type ConstManifold::getView<R, F>(RefMat val, size_t ir, size_t ic) const
  {
    return ptr_->getView<R, F>(val, ir, ic);
  }

  template<>
  inline typename ConstViewReturnType<R, F>::Type ConstManifold::getConstView<R, F>(const ConstRefMat& val, size_t ir, size_t ic) const
  {
    return ptr_->getConstView<R, F>(val, ir, ic);
  }

  template<>
  inline typename ViewReturnType<T, F>::Type ConstManifold::getView<T, F>(RefMat val, size_t ir, size_t ic) const
  {
    return ptr_->getView<T, F>(val, ir, ic);
  }

  template<>
  inline typename ConstViewReturnType<T, F>::Type ConstManifold::getConstView<T, F>(const ConstRefMat& val, size_t ir, size_t ic) const
  {
    return ptr_->getConstView<T, F>(val, ir, ic);
  }

  template<>
  inline typename ViewReturnType<F, F>::Type ConstManifold::getView<F, F>(RefMat val, size_t ir, size_t ic) const
  {
    return ptr_->getView<F, F>(val, ir, ic);
  }

  template<>
  inline typename ConstViewReturnType<F, F>::Type ConstManifold::getConstView<F, F>(const ConstRefMat& val, size_t ir, size_t ic) const
  {
    return ptr_->getConstView<F, F>(val, ir, ic);
  }
}

