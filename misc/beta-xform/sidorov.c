/*
 * sidorov.c
 *
 * Understand bit-sequence mappings.  The expander and compresor
 * functions. This time with alternative expansions, per sidorov.
 *
 * Linas Vepstas Dec 2017; Sept 2020
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// The compressor
double cpr(double y, double K)
{
	// Iterate on y using mashed Bernoulli, and extract symbol dynamics
	char nbits[50];
	for (int i=0; i<50; i++)
	{
		if (0.5 <= y)
			y -= 0.5;
			nbits[i] = 1;
		}
		else nbits[i] = 0;
		y *= 2.0*K;
	}

	// Reconstruct x in a ordinary bernoulli sequence.
	double acc = 1.0e-30;
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

double extcpr(double y, double K)
{
	if (y < K) return cpr(y, K);
	return 2.0 * extcpr(y/(2.0*K), K);
}


// should be same as pdr(cpr(x))
double der(double y, double K)
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

	// Reconstruct x in a mashed bernoulli sequence.
	double acc = 1.0e-30;
	for (int i=0; i<50; i++)
	{
		acc *= 1.0 / (2.0*K);
		if (nbits[50-i-1])
		{
			acc += 0.5;
		}
	}

	return acc;
}

double len(double y, double K)
{
	int longest = 0;
	int len = 0;
	// Iterate on y using mashed Bernoulli, and extract symbol dynamics
	for (int i=0; i<50; i++)
	{
		if (0.5 <= y)
		{
			y -= 0.5;
			len ++;
			if (longest < len) longest = len;
		}
		else
		{
			len = 0;
		}
		y *= 2.0*K;
	}

	return longest;
}

double flen(double K)
{
	return -log(2.0*(1.0-K)) / log (2.0*K);
}

int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K\n", argv[0]);
		exit (1);
	}
	double lam = atof(argv[1]);
	double rng = atof(argv[2]);

#if 1
	int npts = 16803;
	for (int i=0; i<npts; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) npts);
		double y = pdr (x, lam);
		double z = cpr (x, lam);
		printf("%d	%g	%g	%g\n", i, x, y, z);
}
