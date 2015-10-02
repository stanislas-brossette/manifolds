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

#include <manifolds/Manifold.h>
#include <manifolds/defs.h>

namespace mnf
{
  class MANIFOLDS_API ConstSubPoint
  {
  protected:
    ConstSubPoint(ConstManifold M, const ConstRefVec& val);

  public:
    ConstSubPoint(const ConstSubPoint&);

  private:
    ConstSubPoint& operator=(const ConstSubPoint&);

  public:
    ~ConstSubPoint();

    //get the value vector of this point
    ConstRefVec value() const;

    //get a sub point
    const ConstSubPoint operator()(size_t i) const;

    //get the data of a sub point
    //P[i] is equivalent to P(i).value()
    ConstSegment operator[](size_t i) const;

    ConstManifold getManifold() const;
    const Eigen::IOFormat& format() const;

    std::string toString(std::string& prefix, const Eigen::IOFormat& fmt) const; //Dislays point in representation space

  protected:
    ConstManifold manifold_;

    /// \internal We keep value as a non const reference here for easy use in SubPoint. 
    /// This require a const_cast upon building ConstSubPoint, however we honor the
    /// constness in this class.
    RefVec value_;
    mutable Eigen::IOFormat format_;

    friend inline std::ostream& operator<< (std::ostream& os, const ConstSubPoint& x);
  };

  class MANIFOLDS_API SubPoint : public ConstSubPoint
  {
  protected:
    SubPoint(ConstManifold M, RefVec val);

  public:
    SubPoint(const SubPoint&);

  private:
    SubPoint& operator=(const SubPoint&);

  public:
    //get the data the point
    RefVec value();
    ConstRefVec value() const { return ConstSubPoint::value(); }

    //get a sub point
    SubPoint operator()(size_t i);
    ConstSubPoint operator()(size_t i) const { return ConstSubPoint::operator()(i); }

    //get the data of a sub point
    //P[i] is equivalent to P(i).value()
    Segment operator[](size_t i);
    ConstSegment operator[](size_t i) const { return ConstSubPoint::operator[](i); }

  };


  class MANIFOLDS_API PointMemory
  {
  protected:
    PointMemory(Index size);
    PointMemory(const ConstRefVec& v);
    Eigen::VectorXd& getMem();

  private:
    Eigen::VectorXd mem_;
  };


  class MANIFOLDS_API Point : public PointMemory, public SubPoint
  {
    friend class Manifold_Base;
  private:  //only Manifold can create Point
    Point(ConstManifold M);
    Point(ConstManifold M, const ConstRefVec& val);

  public:
    Point(const Point& other);
    Point(const ConstSubPoint& other);

    /// \internal For now, we keep operations on Point only. ConstSubPoint and SubPoint
    /// are only intended for memory read/write, not Manifold operations. 
    Point& increment(const ConstRefVec& v);
    Point& operator=(const Point& x);

    /// \brief Computes a new point that is the result of a retractation of v 
    /// at the current point x. \f$ out = \phi_x(v) \f$
    Point retractation(const ConstRefVec& v) const;
    /// \brief Fills vector out with the value of a point that is the result of a retractation of v 
    /// at the current point x. \f$ out = \phi_x(v) \f$
    void retractation(RefVec out, const ConstRefVec& v) const;
    /// \brief Fills point out with the value of a point that is the result of a retractation of v 
    /// at the current point x. \f$ out = \phi_x(v) \f$
    void retractation(Point& out, const ConstRefVec& v) const;

    /// \brief Computes a vector that is the pseudoLog of point y at the current point x
    /// \f$ out = Log_x(y) \f$
    Eigen::VectorXd pseudoLog(const Point& y) const;
    /// \brief Fills vector out with a vector that is the pseudoLog of point y at the current point x
    /// \f$ out = Log_x(y) \f$
    void pseudoLog(RefVec out, const Point& y) const;

    /// \brief Computes a vector that is the pseudoLog0 of current point x at point 0
    /// \f$ out = Log_0(x) \f$
    Eigen::VectorXd pseudoLog0() const;
    /// \brief Fills a vector as the pseudoLog0 of current point x at point 0
    /// \f$ out = Log_0(x) \f$
    void pseudoLog0(RefVec v) const;
  
    /// \brief returns the typical Magnitude of the manifold
    Eigen::VectorXd typicalMagnitude() const;

    /// \brief returns the Trust Magnitude of the manifold
    Eigen::VectorXd trustMagnitude() const;


    /// \brief Checks that the current point's value represents a valid point of its Manifold
    bool isInM(double prec = 1e-12) const;

    /// \brief Checks that the current point's value represents a valid point of its Manifold
    bool isInTxM(const ConstRefVec& v, const double& prec = 1e-8) const;

    /// \brief Returns the Dimension of its Manifold
    Index getDimM() const;
    /// \brief Returns the Dimension of the tangent space of its Manifold
    Index getTangentDimM() const;
    /// \brief Returns the Dimension of the representation space of its Manifold
    Index getRepresentationDimM() const;

    /// \brief Proxy for Manifold::diffRetractation at current Point
    Eigen::MatrixXd diffRetractation() const;
    /// \brief Proxy for Manifold::applyDiffRetractation at current Point
    void applyDiffRetractation(RefMat out, const ConstRefMat& in) const;
    /// \brief Proxy for Manifold::diffPseudoLog0 at current Point
    Eigen::MatrixXd diffPseudoLog0() const;
    /// \brief Proxy for Manifold::applyDiffPseudoLog0 at current Point
    void applyDiffPseudoLog0(RefMat out, const ConstRefMat& in) const;
    /// \brief Proxy for Manifold::applyTransport at current Point
    void applyTransport(RefMat out, const ConstRefMat& in, const ConstRefVec& v) const;
    /// \brief Proxy for Manifold::applyInvTransport at current Point
    void applyInvTransport(RefMat out, const ConstRefMat& in, const ConstRefVec& v) const;

    virtual const Point& format(const Eigen::IOFormat& fmt) const;
  };

  MANIFOLDS_API Point operator+(const Point& x, const ConstRefVec& v);
  MANIFOLDS_API Eigen::VectorXd operator-(const Point& x, const Point& y);

  inline std::ostream& operator<< (std::ostream& os, const ConstSubPoint& x)
  {
    std::string prefix("");
    os << x.toString(prefix, x.format());
    return os;
  }
}

