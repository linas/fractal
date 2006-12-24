/*
 * mp-gamma.h
 *
 * Compute gamma function for various complex arguments
 *
 * Linas Vepstas December 2006
 */

#include <gmp.h>
#include "mp-complex.h"

/**
 * fp_gamma -- compute Gamma(x)=factorial(x-1) for real argument
 *
 * Uses simple, quickly converging algo-- A&S 6.1.33
 *
 * The caching version skips the calculation, if called again with
 * the same value of ex (up to the nprec precision bits).
 */
void fp_gamma (mpf_t gam, const mpf_t ex, int prec);
void fp_gamma_cache (mpf_t gam, const mpf_t ex, int prec);

/**
 * cpx_gamma -- compute Gamma(x)=factorial(x-1) for complex argument
 *
 * Uses simple, quickly converging algo-- A&S 6.1.33
 *
 * The caching version skips the calculation, if called again with
 * the same value of ex (up to the nprec precision bits).
 */
void cpx_gamma (cpx_t gam, const cpx_t ex, int prec);
void cpx_gamma_cache (cpx_t gam, const cpx_t ex, int prec);
