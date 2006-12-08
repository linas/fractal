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
 */
void fp_gamma (mpf_t gam, const mpf_t ex, int prec);

/**
 * cpx_gamma -- compute Gamma(x)=factorial(x-1) for complex argument
 *
 * Uses simple, quickly converging algo-- A&S 6.1.33
 */
void cpx_gamma (cpx_t gam, const cpx_t ex, int prec);
