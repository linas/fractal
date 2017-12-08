/*
 * remap.c
 *
 * Understand bit-sequence mappings.
 *
 * Linas Vepstas Dec 2017
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double bmap(double x, double Kay)
{
	// Decompose x into a bit sequence.
	char nbits[50];
	for (int i=0; i<50; i++)
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
	double acc = 0.1;
	for (int i=0; i<50; i++)
	{
		acc *= 1.0 / (2.0*Kay);
		if (nbits[50-i-1])
		{
			acc += 0.5;
		}
	}

	return acc;
}

double unbmap(double y, double K)
{
	// Iterate on y using mashed Bernoulli, and extract symbol dynamics
	char nbits[50];
	for (int i=0; i<50; i++)
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
	double acc = 0.1;
	for (int i=0; i<50; i++)
	{
		acc *= 0.5;
		if (nbits[50-i-1])
		{
			acc += 0.5;
		}
	}

	return acc;
}

int main (int argc, char* argv[])
{
	double lam = atof(argv[1]);

	int npts = 803;
	for (int i=0; i<=npts; i++)
	{
		double x = ((double) i)/ ((double) npts);
		double y = unbmap (x, lam);
		double z = bmap (x, lam);
		printf("%d	%g	%g	%g\n", i, x, y, z);
	}
}
