
/*
 * moebius.h
 * 
 * Return the moebius function of an integer.
 * not intended for large integers, works only for small integers
 * due to poor-mans factorization algo.
 *
 * Linas Vepstas Jaunuary 2005
 */

#ifdef   __cplusplus
extern "C" {
#endif

/* classic Moebius mu function */
int moebius_mu (int n);

/* Mertens function, summation of mu */
int mertens_m (int n);

/* The number of prime factors of a number */
int liouville_omega (int n);

/* The Liouville lambda function */
int liouville_lambda (int n);

/* The Mangoldt lambda function */
long double mangoldt_lambda (int n);

#ifdef   __cplusplus
};
#endif
