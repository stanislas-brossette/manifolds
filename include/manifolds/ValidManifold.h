#ifndef _MANIFOLD_VALID_OBJECT_H_
#define _MANIFOLD_VALID_OBJECT_H_

#include <cstdio>
#include <assert.h>

namespace pgs
{
  class ValidManifold
  {
  private:
    static const int FLAG = 1736274519; //a random number
#ifndef NDEBUG
  public:
    ValidManifold() : flag_(FLAG) {}
    ~ValidManifold() { flag_ = -271828182; }

    bool isValid() const { return flag_ == FLAG; }

    bool seeMessageAbove() const 
    {
#  ifndef PGS_ASSERT_THROW
      printf("It appears that you're trying to call a method from a manifold that does not exist \
anymore. Possible cause: the manifold was statically created in a scope (e.g. a function) which was \
left since then. For Stan: https://www.youtube.com/watch?v=Yy8MUnlT9Oo\n"); 
#  endif
      return false; 
    }
  private:
    int flag_;
#endif
  };
}

#endif //_MANIFOLD_VALID_OBJECT_H_
