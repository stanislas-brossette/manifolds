#ifndef _MANIFOLDS_ASSERT_
#define _MANIFOLDS_ASSERT_

#include <cassert>
#include <string>

namespace pgs
{
  struct pgs_exception
  {
  };

#ifdef PGS_ASSERT_THROW 
#  ifndef NDEBUG
  inline void pgs_assert(bool value)
  {
    if(!value)
    {
      throw pgs_exception();
    }
  }
#  else
#    define pgs_assert(expr) assert(expr);
#  endif
#else
#  define pgs_assert(expr) assert(expr);
#endif //_PGS_ASSERT_THROW_
}

#endif //_MANIFOLDS_ASSERT_

