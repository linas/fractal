/* 
 * mp-harmonic.h
 *
 * Compute assorted harmonic functions to arbitrary precision
 *
 * Linas Vepstas October 2006
 */

#include <gmp.h>
#include "mpcomplex.h"

/**
 * mp_pow_rc-- return (k+q)^s for complex s, integer k, real q.
 *
 * If q is held fixed, and k varied, then the values are cached,
 * allowing improved algorithm speeds.
 */

void mp_pow_rc (cpx_t *diri, int k, mpf_t q, cpx_t *ess, int prec);
