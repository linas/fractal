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

/*
 * Exponential generating function for the greatest common factor.
 */
void cpx_gpf_exponential(cpx_t sum, cpx_t z, int prec);

/*
 * Pochhammer-rising generating function for the greatest common factor.
 * This is like the ordinary generating function but uses pochmamer.
 */
void cpx_gpf_poch_rising(cpx_t sum, cpx_t z, int prec);
void cpx_gpf_poch_falling(cpx_t sum, cpx_t z, int prec);

#ifdef  __cplusplus
};
#endif
