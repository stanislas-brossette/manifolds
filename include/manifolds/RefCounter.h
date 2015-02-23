#ifndef _MANIFOLDS_REF_COUNTER_H_
#define _MANIFOLDS_REF_COUNTER_H_

#include <cassert>
#include <manifolds/defs.h>
#include <manifolds/Point.h>
#include <manifolds/pgs_assert.h>

namespace pgs
{
  /// \brief object containing a counter and that cannot be destroyed if the
  /// counter is not at 0.
  class RefCounter
  {
    public:
      RefCounter()
#ifndef NDEBUG
        :count_(0)
#endif
      {
      }

      ~RefCounter()
      {
#ifndef NDEBUG
        pgs_assert(count_ == 0 && "You cannot destroy this manifold because some points still depend on it");
#endif
      }

    protected:
      void incrementRefCounter() const
      {
#ifndef NDEBUG
        count_++;
#endif
      }
      void decrementRefCounter() const
      {
#ifndef NDEBUG
        pgs_assert(count_>0 && "You cannot decrement when no point exist");
        count_--;
#endif
      }

    private:
#ifndef NDEBUG
      mutable int count_;
#endif
      friend void Point::registerPoint();
      friend void Point::unregisterPoint();
  };
}

#endif //_MANIFOLDS_REF_COUNTER_H_

