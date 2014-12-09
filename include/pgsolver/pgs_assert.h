#ifndef _PGS_ASSERT_
#define _PGS_ASSERT_

#include <cassert>
#include <string>

namespace pgs
{
  struct pgs_exception
  {
  };

  inline void pgs_assert(bool value, std::string)
  {
    if(!value)
    {
      //std::cerr << "PGS ASSERT error is " << message << std::endl;
#ifdef PGS_ASSERT_THROW
      throw pgs_exception(); 
#else
      assert(value);
#endif //_PGS_ASSERT_THROW_
    }
  }
    
}

#endif //_PGS_ASSERT_
