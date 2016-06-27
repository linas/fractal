/*
 * FUNCTION:
 * Return greatest prime factor (largest prime divisor)
 *
 * HISTORY:
 * April 2016 -- linas
 * June 2016 -- linas
 */

#ifdef   __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------ */
/* Return the greatest prime factor, 32/64-bit accurate,
 * (depending on size of "unsigned long").
 *
 * Results undefined if input is zero or negative.
 */
unsigned long gpf(unsigned long n);

/* ------------------------------------------------------------ */
/* Return a random number having a distribution similar to the
 * greatest prime factor.  The returned random number is guaranteed
 * to be prime, to be between 2 and n (inclusive) and to occur with
 * roughly the frequency that primes occur in GPF (falling off as
 * 1/k for the k'th prime).  Use srandom() to change the randome seed.
 */
unsigned long pseudo_gpf(unsigned long n);

#ifdef   __cplusplus
};
#endif
