/*
 * isqrt.h
 *
 * Return the integer square-root of an integer.
 *
 * Linas Vepstas October 2016
 */

#ifdef   __cplusplus
extern "C" {
#endif

/**
 * Compute the integer square-root of x.
 * This is an OK algorithm, not as fast as debruijn
 */
int isqrt(int x);
int integer_sqrt(int x);

/**
 * Compute the integer nth-root of x.
 * This is an OK algorithm.
 */
int integer_nth_root(int x, int n);

#ifdef   __cplusplus
};
#endif
