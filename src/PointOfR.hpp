#ifndef _POINT_OF_R_H
#define _POINT_OF_R_H

#include <eigen3/Eigen/Core>
#include "Point.hpp"

namespace pgs
{
  class PointOfR : public Point
  {
    public:
      PointOfR(const int& dim)
        : Point(dim)
      {
      }
    private:
  };
}

#endif
