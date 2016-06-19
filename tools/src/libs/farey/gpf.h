/*
 * FUNCTION:
 * Return greatest prime factor (largest prime divisor)
 *
 * HISTORY:
 * April 2016 -- linas
 */

#ifdef   __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------ */
/*
 * initialize thread-safe locks for the cached version.
 */
void gpf_init(void);

/* Return the greatest prime factor, 32/64-bit accurate,
 * (depending on size of "unsigned long").
 *
 * Results undefined if input is zero.
 */
unsigned long gpf(unsigned long n);

#ifdef   __cplusplus
};
#endif
