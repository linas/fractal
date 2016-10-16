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
 * Ordinary generating function for arithmetic series.
 */
void cpx_ordinary_genfunc(cpx_t sum, cpx_t z, int prec, int (*func)(int));

/*
 * Exponential generating function for arithmetic series.
 */
void cpx_exponential_genfunc(cpx_t sum, cpx_t z, int prec, int (*func)(int));

#ifdef  __cplusplus
};
#endif
