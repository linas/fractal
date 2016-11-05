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
 * Compute p raised to the n'th power. Fast bit-shift algo.
 */

long ipow (long p, long n);

/**
 * Compute the divisor arithmetic function.
 * Returns the number of divisors of n.
 */
long divisor (long int n);

/**
 * Sigma arithmetic series, equals divisor arith series for a=0
 * Computes the divisors of n, raises each to the a'th power, and
 * returns thier sum.
 *
 * Implemented with fast recursive factorization algo.
 *
 * sigmaf is similar, but allows any floating-point exponent.
 * sigmaf is implemented with slow brute force algo.
 */
long sigma (long n, long a);
long double sigmaf (long n, long double a);

/**
 * Much like the sigma arithmetic series, except that an extra
 * log factor is included.   That is, this:
 * Computes the divisors of n, raises each to the a'th power,
 * multiplies the last by logn, and then returns thier sum.
 *
 * sigmalog is implemented with slow brute force algo.
 */
long double sigmalog (long n, long double a);

/**
 * Compute the unitary divisor arithmetic function.
 * A divisor d of is unitary if d and n/d are coprime.
 * Thus, sigma(n,k) = sum_{d|n, gcd(n,n/d)=1} d^k
 * Note: sigma(n,0) = 2^little_omega(n)
 */
long sigma_unitary (long n, long k);

/**
 * Sigma function for k=1; values are cached for performance.
 */
long sigma_one (long n);

/**
 * Parition function
 */
long partition (long n);
unsigned __int128 partitionll (long n);

/** classic Moebius mu function */
long moebius_mu (long n);

/** Mertens function, summatory function of mu */
long mertens_m (long n);

/** Carmichael's lambda function (the funny variant on totient) */
long carmichael_lambda (long n);

/** The number of distinct prime factors of a number. OEIS  A001221 */
long little_omega (long n);

/** The number of prime factors of a number. OEIS  A001222 */
long big_omega (long n);

/** The Liouville lambda function */
long liouville_lambda (long n);

/**
 * The von Mangoldt Lambda function.
 * Returns von Mangoldt Lambda for n, which is
 * log(p) if n is a power of prime p, otherwise
 * returns zero.
 */
long double mangoldt_lambda (long n);
long double mangoldt_lambda_cached (long n);

/** 
 * Exp of the von Mangoldt Lambda function.
 * Returns p if n=p^k for prime p and integer k.
 * Else returns 1.
 */
long  exp_mangoldt_lambda (long n);

/**
 * The indexed von Mangoldt Lambda function
 * Returns the n'th non-zero von Mangoldt value
 */
long double mangoldt_lambda_indexed (long n);
unsigned long mangoldt_lambda_index_point (long n);

/**
 * Thue-Morse
 */
long thue_morse(long n);

#ifdef   __cplusplus
};
#endif
