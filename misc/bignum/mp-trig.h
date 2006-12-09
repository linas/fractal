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
 *
 * The complex exp is built up from the real trig functions.
 * The complex trig functions are built up from the complex exp.
 * In all cases, the most basic idents are used, so these
 * are not speedy!.
 */
void fp_exp (mpf_t ex, const mpf_t z, unsigned int prec);
void fp_sine (mpf_t sine, const mpf_t z, unsigned int prec);
void fp_cosine (mpf_t cosine, const mpf_t z, unsigned int prec);
void cpx_exp (cpx_t ex, const cpx_t z, unsigned int prec);
void cpx_sine (cpx_t sine, const cpx_t z, unsigned int prec);
void cpx_cosine (cpx_t cosine, const cpx_t z, unsigned int prec);
void cpx_tangent (cpx_t tang, const cpx_t z, unsigned int prec);

/**
 * fp_log_m1 -  Floating point logarithm
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */
void fp_log_m1 (mpf_t lg, const mpf_t z, unsigned int prec);
void fp_log (mpf_t lg, const mpf_t z, unsigned int prec);

/**
 * fp_arctan -  Floating point arctangent
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Very slow near y=x
 */
void fp_arctan (mpf_t atn, const mpf_t z, unsigned int prec);
void fp_arctan2 (mpf_t atn, const mpf_t y, const mpf_t x, unsigned int prec);

/*
 * cpx_pow_mpf-- return q^s for complex s, real q.
 *
 * Brute-force algo, this thing is pretty slow, as it requires
 * a logarithm, an exp, sin and cos to be computed, each of which
 * are kinda slow ...
 */
void cpx_pow_mpf (cpx_t powc, const mpf_t q, const cpx_t ess, int prec);

/**
 * fp_pow_rc-- return (k+q)^s for complex s, integer k, real q.
 *
 * If q is held fixed, and k varied, then the values are cached,
 * allowing improved algorithm speeds.
 */

void fp_pow_rc (cpx_t diri, int k, const mpf_t q, const cpx_t ess, int prec);

/* =============================== END OF FILE =========================== */

