/*
 * mp-binomial.h
 *
 * High-precison factorials and binomials, using the 
 * Gnu Multiple-precision library.
 *
 * Also, high-precision values of the series a_n 
 * 
 * Linas Vepstas July 2005
 */

#include <gmp.h>

/* i_poch_rising
 * rising pochhammer symbol, for integer values.
 *
 * Brute force, simple.
 */
void i_poch_rising (mpz_t poch, unsigned int k, unsigned int n);

/** 
 * i_factorial -- the factorial
 */
// #define USE_LOCAL_FACTORIAL
#ifdef USE_LOCAL_FACTORIAL
void i_factorial (mpz_t fact, unsigned int n);
#else
#define i_factorial mpz_fac_ui
#endif /* USE_LOCAL_FACTORIAL */

/* i_binomial
 * Binomial coefficient (n k)
 */
// #define USE_LOCAL_BINOMIAL
#ifdef USE_LOCAL_BINOMIAL
void i_binomial (mpz_t bin, unsigned int n, unsigned int k);
#else 
#define i_binomial mpz_bin_uiui
#endif /* USE_LOCAL_BINOMIAL */

/**
 * i_binomial_sequence -- returns binomial, assumes purely sequential access
 * 
 * This routine assumes that the binomial coefficients will be 
 * accessed in an utterly sequential mode, with k running from 
 * zero to n, and n running from zero to k. For sequential access,
 * this routine is very very fast. Otherwise, random access is used
 * which is considerably slower.
 */
void i_binomial_sequence (mpz_t bin, unsigned int n, unsigned int k);

/** 
 * stirling_first - Stirling Numbers of the First kind, 
 * normalized so that they are all positive.
 * Uses dynamically-sized cache.
 */
void i_stirling_first (mpz_t s, unsigned int n, unsigned int k);
/* A funny off-by-one sum of stirling and binomial */
void i_stirbin_sum (mpz_t s, unsigned int n, unsigned int m);

/* binomial transform of power sum */
void fp_bin_xform_pow (mpf_t bxp, unsigned int n, unsigned int s);

/** 
 * fp_harmonic -- The harmonic number
 */
void fp_harmonic (mpf_t harm, unsigned int n);

/** 
 * fp_poch_rising
 * rising pochhammer symbol (x)_n, for real values of x and integer n.
 *
 * Brute force, simple.
 */
void fp_poch_rising (mpf_t poch, double x, unsigned int n);

/* c_poch_rising
 * rising pochhammer symbol (s)_n, for complex s and integer n.
 *
 * Brute force, simple.
 */
void c_poch_rising_d (mpf_t re_poch, mpf_t im_poch, double re_s, double im_s, unsigned int n);

/* c_poch_rising
 * rising pochhammer symbol (s)_n, for complex s and integer n.
 *
 * Brute force, simple.
 */
void c_poch_rising (mpf_t re_poch, mpf_t im_poch, mpf_t re_s, mpf_t im_s, unsigned int n);

/* fp_binomial
 * Binomial coefficient 
 */
void fp_binomial (mpf_t bin, double s, unsigned int k);

/* c_binomial
 * Complex binomial coefficient
 */
void c_binomial_d (mpf_t re_bin, mpf_t im_bin, double re_s, double im_s, unsigned int k);
void c_binomial (mpf_t re_bin, mpf_t im_bin, mpf_t re_s, mpf_t im_s, unsigned int k);


