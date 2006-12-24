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
 * cpx_polylog -- polylogarithm
 *
 * Li_s(z) = sum_{n=1}^infty z^n/ n^s
 * Works for complex s, z, and uhhhh ,,, this works only for a particular
 * domain that needs to be documented ... may be broken, or incompletely
 * tested.
 */
void cpx_polylog (cpx_t plog, const cpx_t ess, const cpx_t zee, int prec);

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
