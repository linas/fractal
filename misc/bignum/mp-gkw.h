/*
 * mp-gkw.h
 *
 * FUNCTION:
 * Compute matrix elts of the GKW operator.
 *
 * HISTORY:
 * Linas Jan 2010
 */

#include <gmp.h>

#ifdef  __cplusplus
extern "C" {
#endif

// Return the matrix element for the matrix element G_mp of the GKW
// operator, expanded at the x=1 location.
mpf_t gkw(int m, int p);

// Return the continuous-valued version of the GKW operator.
// (the matrix elts occur at integer values)
// This implementation uses GMP multi-precision
mpf_t ache_smooth_mp(double m, double p);

#ifdef  __cplusplus
};
#endif

