/*
 * modular.h
 *
 * FUNCTION:
 * assorted modular functions.
 * -- divisor arithmetic function
 * -- Weierstrass elliptic function g_2 and g_3 invariants
 * -- modular discriminant.
 * -- Klein j-invariant
 *
 * HISTORY:
 * more stuff -- October 2004
 */

extern int modular_max_terms;

/** Compute the divisor arithmetic function
 *  Returns the number of divisors of n.
 */
int divisor (int n);

/** Sigma arithmetic series, equals divisor airth series for a=0 
 *  Computes the divisors of n, raises each to the a'th power, and
 *  returns thier sum.
 */
int sigma (int n, int a);

/* An erdos-borwein-like series sum_n sigma_k(n) x^n/(1-x^n) 
 * Actually computed via a q-series for the divisor function.
 * This is a building block for modular forms.
 */
void erdos_series_c (long double re_q, 
                     long double im_q, 
                     int sa, long double *prep, long double *pimp);

/* Weierstrass elliptic invarient g_2, where q is the nome */
void gee_2_c (long double re_q, long double im_q, 
              long double *pre, long double *pim);

void gee_3_c (long double re_q, long double im_q, 
              long double *pre, long double *pim);

/** The modular discriminant, constructed as g_2^3 - 27g_3^2 . 
 *  Not very accurate, due to the need to take teh difference of two
 *  paritally converged sums.
 */
void disc_c (long double re_q, long double im_q, 
             long double *pre, long double *pim);

/** The modular discriminat, computed as dedekind eta to 24.
 * This version is much better behaved convergent-wise.
 */
void discriminant_c (double re_q, double im_q, double *pre,double *pim);

/* Return euler-product form of the q-series (dedekind eta) */
void euler_prod_c (double re_q, double im_q, double *prep, double *pimp);

/* The dedekind eta multiplies by an additonal factor of q^1/24 */
void dedekind_eta_c (double re_q, double im_q, double *pre, double *pim);

/* The Klein j-invariant g_2^3 / discriminant */
void klein_j_invariant_c (double re_q, double im_q, double *pre, double *pim);

/* --------------------------- END OF FILE ------------------------- */
