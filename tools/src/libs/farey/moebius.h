/*
 * moebius.h
 *
 * Return the Moebius function of an integer.
 *
 * Not intended for large integers, works only for small integers,
 * due to poor-man's factorization algo.
 *
 * Linas Vepstas Jaunuary 2005
 */

#ifdef   __cplusplus
extern "C" {
#endif

/**
 * Compute the integer square-root of x.
 * This is an OK algorithm, not as fast as debruijn
 */
int isqrt(int x);
int integer_sqrt(int x);

/**
 * Compute the integer nth-root of x.
 * This is an OK algorithm.
 */
int integer_nth_root(int x, int n);

/**
 * Compute the divisor arithmetic function.
 * Returns the number of divisors of n.
 */
int divisor (long long int n);

/**
 * Compute the unitary divisor arithmetic function.
 * A divisor d of is unitary if d and n/d are coprime.
 * Thus, this is the sum_{d|n, gcd(n,n/d)=1}
 */
int unitary_divisor (int n);

/** Sigma arithmetic series, equals divisor arith series for a=0
 *  Computes the divisors of n, raises each to the a'th power, and
 *  returns thier sum.
 *  sigmaf is similar, but allows any floating-point exponent.
 */
int sigma (int n, int a);
long double sigmaf (int n, long double a);

/** Much like the sigma arithmetic series, except that an extra
 *  log factor is included.   That is, this:
 *  Computes the divisors of n, raises each to the a'th power,
 *  multiplies the last by logn, and then returns thier sum.
 */
long double sigmalog (int n, long double a);

/** classic Moebius mu function */
int moebius_mu (int n);

/** Mertens function, summatory function of mu */
int mertens_m (int n);

/** Carmichael's lambda function (the funny variant on totient) */
int carmichael_lambda (int n);

/** The number of distinct prime factors of a number. OEIS  A001221 */
int little_omega (int n);

/** The number of prime factors of a number. OEIS  A001222 */
int big_omega (int n);

/** The Liouville lambda function */
int liouville_lambda (int n);

/** The von Mangoldt Lambda function.
 *  Returns von Mangoldt Lambda for n, which is
 *  log(p) if n is a power of prime p, otherwise
 *  returns zero. */
long double mangoldt_lambda (int n);
long double mangoldt_lambda_cached (int n);

/** Exp of the von Mangoldt Lambda function.
 *  Returns p if n=p^k for prime p and integer k.
 *  Else returns 1.
 */
int  exp_mangoldt_lambda (int n);

/** The indexed von Mangoldt Lambda function
 *  Returns the n'th non-zero von Mangoldt value
 *  */
long double mangoldt_lambda_indexed (int n);
unsigned int mangoldt_lambda_index_point (int n);

#ifdef   __cplusplus
};
#endif
