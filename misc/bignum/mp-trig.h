/*
 * mp-trig.h
 *
 * High-precison Elementary functions, using the 
 * Gnu Multiple-precision library.
 *
 * Linas Vepstas July 2005
 */

#include <gmp.h>
#include "mp-complex.h"

/**
 * i_pow - raise n to the m power
 */

void i_pow (mpz_t p, unsigned int n, unsigned int m);

/**
 * fp_inv_pow - raise n to the -m power, where m must be positive. 
 */
void fp_inv_pow (mpf_t p, unsigned int n, unsigned int m);

/**
 * fp_exp -  Floating point exponential
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */
void fp_exp (mpf_t ex, mpf_t z, unsigned int prec);
void fp_sine (mpf_t ex, mpf_t z, unsigned int prec);
void fp_cosine (mpf_t ex, mpf_t z, unsigned int prec);
void cpx_exp (cpx_t *ex, cpx_t *z, unsigned int prec);

/**
 * fp_log_m1 -  Floating point logarithm
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */
void fp_log_m1 (mpf_t lg, mpf_t z, unsigned int prec);
void fp_log (mpf_t lg, mpf_t z, unsigned int prec);

/**
 * fp_arctan -  Floating point arctangent
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */
void fp_arctan (mpf_t atn, mpf_t z, unsigned int prec);

/**
 * fp_pow_rc-- return (k+q)^s for complex s, integer k, real q.
 *
 * If q is held fixed, and k varied, then the values are cached,
 * allowing improved algorithm speeds.
 */

void fp_pow_rc (cpx_t *diri, int k, mpf_t q, cpx_t *ess, int prec);

/* =============================== END OF FILE =========================== */

