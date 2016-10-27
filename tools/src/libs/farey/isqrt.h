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
long isqrt(long x);
long longeger_sqrt(long x);

/**
 * Compute the longeger nth-root of x.
 * This is an OK algorithm.
 */
long longeger_nth_root(long x, long n);

#ifdef   __cplusplus
};
#endif
