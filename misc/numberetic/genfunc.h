/*
 * Generating functions for assorted number-theoretic functions
 * Implementation in bignums.
 *
 * October 2016
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
void cpx_gpf_exponential_shift(cpx_t sum, cpx_t z, int offset, int prec);

#ifdef  __cplusplus
};
#endif
