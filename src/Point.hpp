#ifndef _POINT_H
#define _POINT_H

namespace pgs
{
  class Point
  {
    public:
      Point(const int& dim)
        : dim_(dim)
      {
      }

      int& dim()
      {
        return dim_;
      }
      const int& dim() const
      {
        return dim_;
      }
    private:
      int dim_;
  };
}

#endif
