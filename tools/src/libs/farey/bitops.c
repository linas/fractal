/*
 * bitops.c
 *
 * Bit-wise operations on bit representations of real numbers.
 *
 * Linas Vepstas Dec 2017
 */

#include <math.h>

/**
 * Carry-free multiplication.
 * This "multiplies" two floating point numbers, but then fails
 * to carry the carry bit -- it uses XOR instead of addition.
 */

double mult_xor(double a, double b)
{
	a -= floor(a);
	b -= floor(b);

	// First, create a bit-expansion of a
#define MANTISZ 52
	char abits[MANTISZ];
	for (int i=0; i< MANTISZ; i++)
	{
		if (0.5 <= a) abits[i] = 1;
		else abits[i] = 0;

		a *= 2.0;
		if (1.0 <= a) a -= 1.0;
	}

	// Now, multiply b into a.
	char prod[MANTISZ];
	for (int i=0; i< MANTISZ; i++)
	{
	}
}
