
/** 
 * polylog.h
 *
 * Implement Borwein-style polylogarithm.
 * Also implement the "periodic zeta" and 
 * the Hurwitz zeta function.
 *
 * As of 22 December 2006, seems to be fully functional
 * and correct, and passes tests. The range of convergence
 * is rather limited because of precision/rounding errors.
 *
 * Linas November 2006
 */

#ifndef __POLYLOG_H__
#define __POLYLOG_H__

/**
 * periodic_zeta -- Periodic zeta function 
 *
 * F(s,q) = sum_{n=1}^infty exp(2pi iqn)/ n^s
 *        = Li_s (exp(2pi iq))
 * where 
 * Li_s(z) is the polylogarithm
 *
 * Periodic zeta function is defined as F(s,q) by Tom Apostol, chapter 12
 *
 * Uses the polylog implementation for the calculation,
 * and so is capable of returning results for Re s < 1.
 */
cplex periodic_zeta (cplex s, double q);

/**
 * periodic_zeta_sum -- return value of periodic zeta.
 *
 * Evaluation performed by brute-force sum, thus, 
 * only practical for Re s > 2, otherwise, it will take
 * too long to converge.
 */
cplex periodic_zeta_sum (cplex s, double q);

/**
 * periodic_beta -- Periodic beta function 
 *
 * similar to periodic zeta, but with different normalization
 *
 * beta = 2 Gamma(s+1) (2\pi)^{-s} F(s,q)
 *
 * As of 22 December, seems to be passing the tests -- 
 * that is, it gives the Bernoulli polynomials for integer s,
 * with all the right scale factors and signs, etc. Yay!
 */
cplex periodic_beta (cplex s, double q);

/**
 * hurwitz_zeta -- Hurwitz zeta function
 *
 * Built up from the periodic beta
 * As of 22 December 2006, seems to be working, with some trouble
 * due to insufficient precision.
 */
cplex hurwitz_zeta (cplex ess, double q);

#endif /* __POLYLOG_H__ */
