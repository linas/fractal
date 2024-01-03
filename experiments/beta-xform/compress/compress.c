/*
 * compress.c
 *
 * Common implementation for compressor, expander functions.
 *
 * Linas Vepstas Jan 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// The expander. A binary bit sequence for x is obtained, and
// resummed using K.  Recall K=beta/2.
double pdr(double x, double Kay)
{
#define NBITS 50
	// Decompose x into a bit sequence.
	char nbits[NBITS];
	for (int i=0; i<NBITS; i++)
	{
		if (0.5 <= x)
		{
			x -= 0.5;
			nbits[i] = 1;
		}
		else nbits[i] = 0;
		x *= 2.0;
	}

	// Reconstruct x in a mashed bernoulli sequence.
	double acc = 1.0e-30;
	for (int i=0; i<NBITS; i++)
	{
		acc *= 1.0 / (2.0*Kay);
		if (nbits[NBITS-i-1])
			acc += 0.5;
	}

	return acc;
}

// The compressor. A beta-expansion for y is obtained, i.e. a sequence
// of bits in base-beta, and then resassembled in base two.
// Recall K = beta/2.
double cpr(double y, double K)
{
	// Iterate on y using mashed Bernoulli, and extract symbol dynamics
	char nbits[NBITS];
	for (int i=0; i<NBITS; i++)
	{
		if (0.5 <= y)
		{
			y -= 0.5;
			nbits[i] = 1;
		}
		else nbits[i] = 0;
		y *= 2.0*K;
	}

	// Reconstruct x in a ordinary bernoulli sequence.
	double acc = 1.0e-30;
	for (int i=0; i<NBITS; i++)
	{
		acc *= 0.5;
		if (nbits[NBITS-i-1])
			acc += 0.5;
	}

	return acc;
}
