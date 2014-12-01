#ifndef _PGS_REF_COUNTER_H_
#define _PGS_REF_COUNTER_H_

#include <cassert>

#define DEBUG

namespace pgs
{
  class RefCounter
  {
    public:
      RefCounter()
#ifdef DEBUG
        :count(0)
#endif
      {
      }

      ~RefCounter()
      {
#ifdef DEBUG
        assert(count == 0 && "You cannot destroy this manifold because some points still depend on it");
#endif
      }

    protected:
      void incrementRefCounter() const
      {
#ifdef DEBUG
        count++;
#endif
      }
      void decrementRefCounter() const
      {
#ifdef DEBUG
        assert(count>0 && "You cannot decrement when no point exist");
        count--;
#endif
      }

    private:
#ifdef DEBUG
      mutable int count;
#endif
      friend Point::Point(const Point&);
      friend Point::~Point();
  };
}

#endif //_PGS_REF_COUNTER_H_
