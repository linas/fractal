/*
 * mp_zeta.h
 *
 * High-precison Riemann zeta function, using the 
 * Gnu Multiple-precision library.
 *
 * Also, high-precision values of the series a_n 
 * 
 * Linas Vepstas July 2005
 */

#include <gmp.h>
#include "mp-complex.h"

/* Fixed-point bernoulli number */
void q_bernoulli (mpq_t bern, int n);

/* Compute and return the "exact" result for the zeta function for 
 * any value of even n. Computed top prec decimal places.
 * works by computing the Bernoulli number first.
 */
void fp_zeta_even (mpf_t zeta, unsigned int n, int prec);

/* fp_zeta
 * Floating-point-valued Riemann zeta for positive integer arguments 
 * return value placed in the arg "zeta".
 *
 * Carries out math to prec decimal digits
 */
void fp_zeta (mpf_t zeta, unsigned int s, int prec);

/* Same, using Helmut Hasse convergent algo */
void fp_hasse_zeta (mpf_t zeta, unsigned int s, int prec);

/* Same, using P. Borwein convergent algo */
void fp_borwein_zeta (mpf_t zeta, unsigned int s, int prec);

/* P. Borwein zeta for complex argument */
void fp_borwein_zeta_c (cpx_t zeta, cpx_t ess, int prec);

/* Brute-force summation */
void fp_zeta_brute (mpf_t zeta, unsigned int s, int prec);

/* Stieltjes constants */
// void stieltjes_gamma (mpf_t gam, int n);

/* 
 * Compute a_sub_n
 * the w argument is for the power bit -- 
 */
void a_sub_n (mpf_t a_n, mpf_t w, unsigned int n, unsigned int prec);
void b_sub_n (mpf_t b_n, unsigned int n, unsigned int prec);

/* compute a_sub_s for complex-valued s
 */
void a_sub_s (mpf_t re_a, mpf_t im_a, double re_s, double im_s, unsigned int prec);
void b_sub_s_d (mpf_t re_a, mpf_t im_a, double re_s, double im_s, 
              unsigned int prec, int nterms, double eps);
void b_sub_s (mpf_t re_a, mpf_t im_a, mpf_t re_s, mpf_t im_s, 
              unsigned int prec, int nterms, double eps);

