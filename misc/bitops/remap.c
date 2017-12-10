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

// The expander
double pdr(double x, double Kay)
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

// the compressor
double cpr(double y, double K)
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

double eig(double y, double K)
{
	// if (1.0 < y) return eig(y-1.0, K);
	// if (0.99*K < y && y < 1.01*K) return -0.5;

	if (1.0 < y) {
		// Attempted unsuccessful patch.
		return eig((y-1.0)/K+K, K) ; // - eig(y/(2.0*K), K);
	}
	return cpr(y,K) - 0.5;
}

int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K\n", argv[0]);
		exit (1);
	}
	double lam = atof(argv[1]);

	int npts = 803;
	for (int i=0; i<=npts; i++)
	{
		double x = ((double) i)/ ((double) npts);
#ifdef BASICS
		double y = pdr (x, lam);
		double z = cpr (x, lam);
		printf("%d	%g	%g	%g\n", i, x, y, z);
#endif
// #define SELF_SIM
#ifdef SELF_SIM
		double y = pdr (x, lam);
		double lo = pdr (0.5*x, lam);
		double hi = pdr (0.5+0.5*x, lam);
		double z = pdr (1-x, lam);
		printf("%d	%g	%g	%g	%g	%g\n", i, x, y, lo, hi, z);
#endif
// #define CPR_SIM 1
#if CPR_SIM
		double y = cpr (x, lam);
		double lo = cpr (0.5*x, lam);
		double hi = cpr (0.5+0.5*x, lam);
		double z = cpr (lam*x, lam);
		printf("%d	%g	%g	%g	%g	%g\n", i, x, y, lo, hi, z);
#endif
		double y = eig (x, lam);
		double lo = eig (x/(2.0*lam), lam);
		double hi = eig ((lam+x)/(2.0*lam), lam);
		printf("%d	%g	%g	%g	%g\n", i, x, y, lo, hi);
	}
}
