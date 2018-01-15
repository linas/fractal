/*
 * bitops.c
 *
 * Bit-wise operations on bit representations of real numbers.
 *
 * Linas Vepstas Dec 2017
 */

#include <math.h>

#include "bitops.h"

#define MANTISZ 58

// Convert float to bit string
static void float_to_bits(double a, char* bits)
{
	for (int i=0; i< MANTISZ; i++)
	{
		if (0.5 <= a) bits[i] = 1;
		else bits[i] = 0;

		// shift left one bit
		a *= 2.0;
		if (1.0 <= a) a -= 1.0;
	}
}

// Convert bits to double.
static double bits_to_float(char* bits)
{
	double p = 0.0;
	double h = 0.5;
	for (int i=0; i< MANTISZ; i++)
	{
		if (0 != bits[i]) p += h;
		h *= 0.5;
	}
	return p;
}

/**
 * Carry-free addition.
 * This takes the XOR of the binary representation of two floating
 * point numbers, and returns the result as a float.  Basically, its
 * just ordinary addition, with a failure to to carry the carry bit.
 */

double add_xor(double a, double b)
{
	a -= floor(a);
	b -= floor(b);

	// XXX This is the simple sloppy implementation, not the
	// high-performance one!
	// First, create a bit-expansion of a
	char abits[MANTISZ];
	float_to_bits(a, abits);

	char bbits[MANTISZ];
	float_to_bits(b, bbits);

	// Take the xor
	for (int i=0; i< MANTISZ; i++)
	{
		abits[i] = abits[i] ^ bbits[i];
	}

	// Now convert back to double.
	return bits_to_float(abits);
}

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
	char abits[MANTISZ];
	float_to_bits(a, abits);

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
			// Accumulate (add) a/2^i into the product
			for (int j=i; j<MANTISZ; j++)
			{
				prod[j] = prod[j] ^ abits[j-i];
			}
		}

		// Shift left one bit
		b *= 2.0;
		if (1.0 <= b) b -= 1.0;
	}

	// Now convert back to double.
	return bits_to_float(prod);
}

double multiply(double a, double b)
{
	a -= floor(a);
	b -= floor(b);

	// First, create a bit-expansion of a
	char abits[MANTISZ];
	float_to_bits(a, abits);

	// The sum of the carry bits
	char sum[MANTISZ];
	for (int i=0; i< MANTISZ; i++)
	{
		sum[i] = 0;
	}

	// Now, multiply b into a.
	for (int i=0; i< MANTISZ; i++)
	{
		if (0.5 <= b)
		{
			// Accumulate a/2^i into the product
			for (int j=i; j<MANTISZ; j++)
			{
				if (abits[j-i]) sum[j] ++;
			}
		}

		// Shift left one bit
		b *= 2.0;
		if (1.0 <= b) b -= 1.0;
	}

	// The carry bits 
	char carry[MANTISZ+2];
	carry[MANTISZ+1] = 0;
	for (int i= MANTISZ; 0 < i; i--)
	{
		carry[i] = sum[i-1] + (carry[i+1] >> 1);
	}
	carry[0] = carry[1] >> 1;

	// The product 
	char prod[MANTISZ];
	for (int i=0; i< MANTISZ; i++)
	{
		prod[i] = carry[i] % 2;
	}

	// Now convert back to double.
	return bits_to_float(prod);
}

// #define UNIT_TEST
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

	printf ("its %g otimes %g = %g\n", x,y, mult_xor(x,y));
	double mxy = multiply (x,y);
	printf ("%g x %g = %g = %g delta=%g\n", x,y, mxy, x*y, mxy-x*y);
}
#endif
