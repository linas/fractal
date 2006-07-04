
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

/** classic Moebius mu function */
int moebius_mu (int n);

/** Mertens function, summatory function of mu */
int mertens_m (int n);

/** The number of prime factors of a number */
int liouville_omega (int n);

/** The Liouville lambda function */
int liouville_lambda (int n);

/** The von Mangoldt Lambda function 
 *  Returns von Mangoldt Lambda for n, which is
 *  log(p) if n is a power of prime p, otherwise 
 *  returns zero. */
long double mangoldt_lambda (int n);
long double mangoldt_lambda_cached (int n);

/** The indexed von Mangoldt Lambda function 
 *  Returns the n'th non-zero von Mangoldt value
 *  */
long double mangoldt_lambda_indexed (int n);
unsigned int mangoldt_lambda_index_point (int n);

#ifdef   __cplusplus
};
#endif
