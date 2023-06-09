/*
 * FUNCTION:
 * Return greatest common factor (greatest common divisor)
 *
 * HISTORY:
 * Oct 2004 -- linas
 */

#ifdef   __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------ */
/**
 * Return the greatest common factor, 64-bit accurate,
 * (depending on size of "unsigned long").  AKA gcd.
 *
 * Results undefined if input is zero.
 */
unsigned long gcf64 (unsigned long nume, unsigned long denom);

/**
 * Return the least common multiple, 64-bit accurate,
 * (depending on size of "unsigned long").
 *
 * Note that lcm(a,b) = ab / gcf(a,b)
 */
unsigned long lcm64 (unsigned long a, unsigned long b);

#ifdef   __cplusplus
};
#endif
