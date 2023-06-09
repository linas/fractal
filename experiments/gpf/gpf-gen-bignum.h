/*
 * Generating functions for greatest prime factors.
 * Implementation in bignums.
 *
 * April 2016
 */

#include <mp-complex.h>

#ifdef  __cplusplus
extern "C" {
#endif

/*
 * Ordinary generating function for the greatest common factor.
 */
void cpx_gpf_ordinary(cpx_t sum, cpx_t z, int prec);
void cpx_gpf_ordinary_shift(cpx_t sum, cpx_t z, int offset, int prec);
void cpx_gpf_ordinary_derivative(cpx_t sum, cpx_t z, int offset, int prec);
void cpx_gpf_ordinary_recip(cpx_t sum, cpx_t z, int prec);
// void cpx_gpf_ordinary_s(cpx_t sum, cpx_t z, cpx_t s, int prec);
// void cpx_gpf_ordinary_s(cpx_t sum, cpx_t z, int s, int prec);

/*
 * Exponential generating function for the greatest common factor.
 */
void cpx_gpf_exponential(cpx_t sum, cpx_t z, int prec);
void cpx_gpf_exponential_shift(cpx_t sum, cpx_t z, int offset, int prec);
void cpx_gpf_exponential_newton(cpx_t sum, cpx_t z, int offset, int prec);
void cpx_gpf_exponential_recip(cpx_t sum, cpx_t z, int prec);
void cpx_gpf_exponential_s(cpx_t sum, cpx_t z, cpx_t s, int prec);

/*
 * Like the above, but only the odd/even terms, and an alternating sign.
 */
void cpx_gpf_sine(cpx_t sum, cpx_t z, int prec);
void cpx_gpf_cosine(cpx_t sum, cpx_t z, int prec);
/*
 * Dirichlet generating function for the greatest common factor.
 */
void cpx_gpf_dirichlet(cpx_t sum, cpx_t s, int prec);

/*
 * Pochhammer-rising generating function for the greatest common factor.
 * This is like the exponential generating function, but uses pochhammer
 * instead of x^n. Actually, uses binomial coefficient ... two
 * factorials are needed for convergence.
 */
void cpx_gpf_poch_rising(cpx_t sum, cpx_t z, int prec);
void cpx_gpf_poch_falling(cpx_t sum, cpx_t z, int prec);

void cpx_random_exponential_shift(cpx_t sum, cpx_t z, int offset, int prec);
#ifdef  __cplusplus
};
#endif
