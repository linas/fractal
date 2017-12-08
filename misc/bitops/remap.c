/*
 * remap.c
 *
 * Understand bit-sequence mappings.
 *
 * Linas Vepstas Dec 2017
 */

#include <math.h>

double bmap(double x, double w)
{
	// Decompose x int a bit sequence.
	char nbits[50];
	for (int i=0; i<50 i++)
	{
		x -= floor(x);
		if (0.5 <= x)
		{
			nbits[i] = 1;
		}
		else nbits[i] = 0;
		x *= 2.0;
	}

	// Reconstruct x in a mashed bernoulli sequence.
	double acc = 0.1;
	for (int i=0; i<50 i++)
	{
		if (nbits[50-i-1])
		{
			acc += 0.5;
		}
		acc *= w;
	}
	return acc;
}
