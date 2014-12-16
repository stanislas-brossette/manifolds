#ifndef _LEX_ASSERT_H_
#define _LEX_ASSERT_H_


struct lex_assertion_fired {};

#undef lex_assert

#ifndef NDEBUG
# ifndef LEX_ASSERT_ACTIVE
#  define LEX_ASSERT_ACTIVE
# endif
#endif

#ifdef LEX_ASSERT_ACTIVE

# include <cstdio>
# define lex_assert(xpr) if(!(xpr)) { fprintf(stderr, #xpr"\n"); fflush(stderr); throw lex_assertion_fired(); }

#else // LEX_ASSERT_ACTIVE not defined

# define lex_assert(xpr)     ((void)0)

#endif // LEX_ASSERT_ACTIVE

#endif // _LEX_ASSERT_H_