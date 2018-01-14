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
 * Take the XOR of the binary representatation of two floats.
 * Return the result, converted back into a float.  This is a lot
 * like addition, but failing to carry the carry bit.
 */
double add_xor(double a, double b);

/**
 * Carry-free multiplication.
 * This "multiplies" two floating point numbers, but then fails
 * to carry the carry bit -- it uses XOR instead of addition.
 */
double mult_xor(double a, double b);

#ifdef   __cplusplus
};
#endif
