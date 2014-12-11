#ifndef _PGS_REF_COUNTER_H_
#define _PGS_REF_COUNTER_H_

#include <cassert>
#include <pgsolver/Point.h>

#define DEBUG

namespace pgs
{
  class RefCounter
  {
    public:
      RefCounter()
#ifdef DEBUG
        :count_(0)
#endif
      {
      }

      ~RefCounter()
      {
#ifdef DEBUG
        assert(count_ == 0 && "You cannot destroy this manifold because some points still depend on it");
#endif
      }

    protected:
      void incrementRefCounter() const
      {
#ifdef DEBUG
        count_++;
#endif
      }
      void decrementRefCounter() const
      {
#ifdef DEBUG
        assert(count_>0 && "You cannot decrement when no point exist");
        count_--;
#endif
      }

    private:
#ifdef DEBUG
      mutable int count_;
#endif
      friend void Point::registerPoint();
      friend void Point::unregisterPoint();
  };
}

#endif //_PGS_REF_COUNTER_H_
