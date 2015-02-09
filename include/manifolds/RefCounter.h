#ifndef _MANIFOLDS_REF_COUNTER_H_
#define _MANIFOLDS_REF_COUNTER_H_

#include <cassert>
#include <manifolds/defs.h>
#include <manifolds/Point.h>

#define DEBUG

namespace pgs
{
  /// \brief object containing a counter and that cannot be destroyed if the
  /// counter is not at 0.
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

#endif //_MANIFOLDS_REF_COUNTER_H_

