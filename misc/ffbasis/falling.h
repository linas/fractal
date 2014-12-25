/**
 * falling.h
 *
 * Falling factorial basis for the divisor function
 *
 * Linas Vepstas December 2014
 */

/**
 * The operator that transforms binomial (falling factorial) into a
 * Dirichlet series. Defined as:
 *    E_mk = (1/m - 1)^k
 *
 * We assume m>0 and that k >= 0
 */
long double E_mk(unsigned int m, unsigned int k);

/**
 * A vector in the kernel of E aka a_k
 * If we did this correctly, a_k is given by solving
 *    sum_k=0^infty a_k x^k = sin (2pi/(1+x))
 * This returns a_k.
 */
long double a_k(unsigned int k);
