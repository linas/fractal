/** 
 * mp-polylog.h
 *
 * Implement Borwein-style polylogarithm.
 * Also implement the "periodic zeta" and 
 * the Hurwitz zeta function.
 *
 * Linas November 2006
 */
#include "mp-complex.h"

/**
 * cpx_polylog_nint -- compute the polylogarithm at negetive integers
 *
 * Li_{-n}(z) 
 * At the negative integers, the polylog is a rational function,
 * meromorphic everywhere except for multiple poles at z=1.
 */
void cpx_polylog_nint (cpx_t plog, unsigned int negn, const cpx_t zee);

/**
 * cpx_polylog_sum -- compute the polylogarithm by direct summation
 *
 * Li_s(z) = sum_{n=1}^infty z^n/ n^s
 * 
 * The magnitude of z must be less than one in order for the 
 * summation to be caqrried out.
 *
 * Caches intermediate results, so that overall performance is
 * considerably better if z is varied while s is held fixed.
 */
void cpx_polylog_sum (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec);

/**
 * cpx_polylog -- polylogarithm
 *
 * Li_s(z) = sum_{n=1}^infty z^n/ n^s
 * 
 * Works for general complex s, z; lightly tested, may be buggy.
 * Watch out for branchpoint at z=1.
 *
 * Returns a non-zero value if algo was unable to evaluate at
 * the given point.
 */
int cpx_polylog (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec);

/**
 * cpx_periodic_zeta -- Periodic zeta function 
 *
 * F(s,q) = sum_{n=1}^infty exp(2pi iqn)/ n^s
 *        = Li_s (exp(2pi iq))
 * where 
 * Li_s(z) is the polylogarithm
 *
 * Periodic zeta function is defined as F(s,q) by Tom Apostol, chapter 12
 */
void cpx_periodic_zeta (cpx_t z, const cpx_t ess, const mpf_t que, int prec);

/**
 * cpx_periodic_beta -- Periodic beta function 
 *
 * Similar to periodic zeta, but with different normalization
 *
 * beta = 2 Gamma(s+1) (2\pi)^{-s} F(s,q)
 *
 * Caches intermediate terms, and so performance is much better 
 * if s is held const, while q is varied.
 */
void cpx_periodic_beta (cpx_t zee, const cpx_t ess, const mpf_t que, int prec);

/**
 * cpx_hurwitz_zeta -- Hurwitz zeta function
 * Returns zeta = sum_{n=0}^infty 1/(n+q)^s
 * Accepts complex s, real-valued q.
 *
 * Built up from the fast polylogarithm algo
 * Caches intermediate terms, and so performance is much better 
 * if s is held const, while q is varied.
 */
void cpx_hurwitz_zeta (cpx_t hzeta, const cpx_t ess, const mpf_t que, int prec);
