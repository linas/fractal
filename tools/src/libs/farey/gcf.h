
/* 
 * FUNCTION:
 * return greatest common factor (greatest common divisor)
 *
 * HISTORY:
 * Oct 2004 -- linas
 */

#ifdef   __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------ */
/* Return the greatest common factor, 32-bit in accurate.
 * Assumes that inputs are positive integers; results undefined 
 * for neg integers or zero.
 */
int gcf32 (int nume, int denom);

#ifdef   __cplusplus
};
#endif
