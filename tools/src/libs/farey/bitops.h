/*
 * bitops.h
 *
 * Bit-wise operations on bit representations of real numbers.
 *
 * Linas Vepstas Dec 2017
 */

#ifdef   __cplusplus
extern "C" {
#endif

/**
 * Carry-free multiplication.
 * This "multiplies" two floating point numbers, but then fails
 * to carry the carry bit -- it uses XOR instead of addition.
 */
double mult_xor(double a, double b);

#ifdef   __cplusplus
};
#endif
