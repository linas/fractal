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

void fp_prt (char * str, mpf_t val);

/* i_poch_rising
 * rising pochhammer symbol, for integer values.
 *
 * Brute force, simple.
 */
void i_poch_rising (mpz_t poch, unsigned int k, unsigned int n);

/* i_factorial -- the factorial
 */
void i_factorial (mpz_t fact, unsigned int n);

/* i_binomial
 * Binomial coefficient (n k)
 */
void i_binomial (mpz_t bin, unsigned int n, unsigned int k);

/* fp_poch_rising
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

void c_poch_rising (mpf_t re_poch, mpf_t im_poch, double re_s, double im_s, unsigned int n);

/* fp_binomial
 * Binomial coefficient 
 */
void fp_binomial (mpf_t bin, double s, unsigned int k);

/* c_binomial
 * Complex binomial coefficient
 */
void c_binomial (mpf_t re_bin, mpf_t im_bin, double re_s, double im_s, unsigned int k);

/* fixed-point bernoulli number */
void q_bernoulli (mpq_t bern, int n);

/* Floating point exponential */
void fp_exp (mpf_t ex, mpf_t z, unsigned int prec);

/* fp_euler
 * return Euler-Mascheroni const
 */
void fp_euler_mascheroni (mpf_t gam);

void fp_pi (mpf_t pi);

/* return e^pi */
void fp_e_pi (mpf_t e_pi);

void fp_zeta2 (mpf_t zeta);
void fp_zeta3 (mpf_t zeta);
void fp_zeta5 (mpf_t zeta);
void fp_zeta7 (mpf_t zeta);
void fp_zeta9 (mpf_t zeta);

/* Compute and return the "exact" result for the zeta function for 
 * any value of even n 
 */
void fp_zeta_even (mpf_t zeta, unsigned int n);



/* fp_zeta
 * Floating-point-valued Riemann zeta for positive integer arguments 
 * return value placed in the arg "zeta".
 *
 * Simple-minded algo, carries out math to prec decimal digits
 *
 * OLEIS A046988
 *  Numerator of zeta(2n)/Pi^(2n)
 * 0,1,1,1,1,1,691,2,3617,43867,174611,155366,236364091,
 * 1315862,6785560294,6892673020804,7709321041217,151628697551,
 * 26315271553053477373,308420411983322,261082718496449122051,
 * 3040195287836141605382,5060594468963822588186
 *
 * OLEIS A002432 
 *  Denominator of zeta(2n)/Pi^(2n)
 *  6,90,945,9450,93555,638512875,18243225,325641566250,
 *  38979295480125,1531329465290625,13447856940643125,
 *  201919571963756521875,11094481976030578125,
 *  564653660170076273671875,5660878804669082674070015625,
 *  62490220571022341207266406250,12130454581433748587292890625
 */
void fp_zeta (mpf_t zeta, unsigned int s, int prec);

/* 
 * Compute a_sub_n
 * the w argument is for the power bit -- 
 */
void a_sub_n (mpf_t a_n, mpf_t w, unsigned int n, unsigned int prec);

/* compute a_sub_s for complex-valued s
 */
void a_sub_s (mpf_t re_a, mpf_t im_a, double re_s, double im_s, unsigned int prec);

