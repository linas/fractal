/*
 * dirichlet.h
 *
 * Return iterated dirichlet convolution.
 *
 * Linas Vepstas January 2019
 */

#ifdef   __cplusplus
extern "C" {
#endif

/**
 * Return iterated Dirichlet convolution of the unit function,
 * and it's inverse (the Moebius functions).
 * unit(k,n) == k'th iteration, n'th value. (for positive n, only)
 * k may be any integer (pos, neg, zero)
 *
 * So:
 * -1 == k is the mobius mu(n)  -- unit(-1, n) == mobius_mu(n)
 * +1 == k is all-ones          -- unit(1, n) == 1 for all n
 *  0 == k is the identity:     -- unit(0, 1)=1 otherwise unit(0, n)=0 for n>1
 *
 * Otherwise its recursively defined:
 * k in general, unit(k) = unit(k-d) * unit(d) for all integers d
 * where * is the Dirichlet convolution.
 */
long unit(long k, long n);

#ifdef   __cplusplus
};
#endif
