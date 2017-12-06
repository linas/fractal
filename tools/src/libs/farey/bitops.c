/*
 * bitops.c
 *
 * Bit-wise operations on bit representations of real numbers.
 *
 * Linas Vepstas Dec 2017
 */

#include <math.h>

#include "bitops.h"

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

		// shift left one bit
		a *= 2.0;
		if (1.0 <= a) a -= 1.0;
	}

	// The carry-free product
	char prod[MANTISZ];
	for (int i=0; i< MANTISZ; i++)
	{
		prod[i] = 0;
	}

	// Now, multiply b into a.
	for (int i=0; i< MANTISZ; i++)
	{
		if (0.5 <= b)
		{
			// accumulate (add) a/2^i into the product
			for (int j=i; j<MANTISZ; j++)
			{
				prod[j] = prod[j] ^ abits[j-i];
			}
		}

		// shift left one bit
		b *= 2.0;
		if (1.0 <= b) b -= 1.0;
	}

	// Now convert back to double.
	double p = 0.0;
	double h = 0.25;
	for (int i=0; i< MANTISZ; i++)
	{
		if (0 != prod[i]) p += h;
		h *= 0.5;
	}

	return p;
}

#ifdef UNIT_TEST

#include <stdio.h>
#include <stdlib.h>

int main (int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "%s: expecting two float arguments\n", argv[0]);
		exit(1);
	}
	double x = atof(argv[1]);
	double y = atof(argv[2]);

	printf ("its %g x %g = %g\n", x,y, mult_xor(x,y));
}
#endif
