/*
 * Generating functions for greatest prime factors.
 * Implementation in bignums.
 *
 * April 2016
 */

#include <mp-complex.h>

/*
 * Ordinary generating function for the greatest common factor.
 */
void cpx_gpf_ordinary(cpx_t sum, cpx_t z, int prec);

/*
 * Exponential generating function for the greatest common factor.
 */
void cpx_gpf_exponential(cpx_t sum, cpx_t z, int prec);
