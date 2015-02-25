#ifndef _MANIFOLDS_POINT_H_
#define _MANIFOLDS_POINT_H_

#include <manifolds/defs.h>

namespace pgs
{
  class Manifold;

  class MANIFOLDS_API ConstSubPoint
  {
  protected:
    ConstSubPoint(const Manifold& M, RefVec val);
    ConstSubPoint(const Manifold& M, const ConstRefVec& val);

  public:
    ConstSubPoint(const ConstSubPoint&);

  private:
    ConstSubPoint& operator=(const ConstSubPoint&);

  public:
    ~ConstSubPoint();

    //get the data of this point
    ConstRefVec value() const;

    //get a sub point
    const ConstSubPoint operator()(size_t i) const;

    //get the data of a sub point
    //P[i] is equivalent to P(i).value()
    ConstSegment operator[](size_t i) const;

    const Manifold& getManifold() const;

    std::string toString(std::string& prefix) const; //Dislays point in representation space

  private:
    void registerPoint();
    void unregisterPoint();

  protected:
    const Manifold& manifold_;

    /// \internal We keep value as a non const reference here for easy use in SubPoint. 
    /// This require a const_cast upon building ConstSubPoint, however we honor the
    /// constness in this class.
    RefVec value_;

    friend inline std::ostream& operator<< (std::ostream& os, const ConstSubPoint& x);
    friend class RefCounter;
  };

  class MANIFOLDS_API SubPoint : public ConstSubPoint
  {
  protected:
    SubPoint(const Manifold& M, RefVec val);

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
    //using ConstSubPoint::operator();

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
  private:  //only Manifold can create Point
    Point(const Manifold& M);
    Point(const Manifold& M, const ConstRefVec& val);

  public:
    Point(const Point& other);
    Point(const ConstSubPoint& other);

    /// \internal For now, we keep operations on Point only. ConstSubPoint and SubPoint
    /// are only intended for memory read/write, not Manifold operations. 
    Point& increment(const ConstRefVec& v);
    Point & operator=(const Point& x);

    friend class Manifold;
  };

  MANIFOLDS_API Point operator+(const Point& x, const ConstRefVec& v);
  MANIFOLDS_API Eigen::VectorXd operator-(const Point& x, const Point& y);

  inline std::ostream& operator<< (std::ostream& os, const ConstSubPoint& x)
  {
    std::string prefix("");
    os << x.toString(prefix);
    return os;
  }
}

#endif //_MANIFOLDS_POINT_H_

