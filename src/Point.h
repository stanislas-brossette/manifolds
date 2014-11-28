#ifndef _PGS_POINT_H_
#define _PGS_POINT_H_

#include "defs.h"

namespace pgs
{
  class Manifold;

  class Point
  {
  private:  //only Manifold can create Point
    Point(const Manifold& M);
    Point(const Manifold& M, const Eigen::VectorXd& val);

  public:
    Point& increment(const Eigen::VectorXd& v);

    //get a sub point
    Point operator()(size_t i) const;

    //get the data of a sub point
    const Eigen::VectorXd& value() const;
    ConstSegment operator[](size_t i) const;
    Segment operator[](size_t i);

    const Manifold& getManifold() const;
    
    std::string toString(std::string& prefix) const; //Dislays point in representation space

    friend inline std::ostream& operator<< (std::ostream& os, const Point& x);
    
  private:
    

  private:
    const Manifold& manifold_;
    Eigen::VectorXd value_;

    friend class Manifold;
  };

  Point operator+(const Point& x, const Eigen::Ref<const Eigen::VectorXd>& v);

  inline std::ostream& operator<< (std::ostream& os, const Point& x)
  {
    std::string prefix("");
    os << x.toString(prefix);
    return os;
  }
}

#endif //_PGS_POINT_H_

