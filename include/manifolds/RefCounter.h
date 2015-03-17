#ifndef _MANIFOLDS_REF_COUNTER_H_
#define _MANIFOLDS_REF_COUNTER_H_

#include <manifolds/defs.h>
#include <manifolds/Point.h>
#include <manifolds/pgs_assert.h>

//180021114 is a version number for which noexcept works, not sure about noexcept(false)
#if defined(_MSC_FULL_VER) && _MSC_FULL_VER < 180021114
#define NOEXCEPT(x)
#else
#define NOEXCEPT(x) noexcept(x)
#endif

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

      // In C++ 11, by default, a destructor cannot launch an exception, setting noexcet to false allows it.
      ~RefCounter() NOEXCEPT(false)
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
      friend void ConstSubPoint::registerPoint();
      friend void ConstSubPoint::unregisterPoint();
  };
}

#endif //_MANIFOLDS_REF_COUNTER_H_

