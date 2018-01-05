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
	double acc = 1.0e-30;
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


// deconstruct/reconstruct
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

// the tent compressor
double tent_cpr(double y, double K)
{
	// Iterate on y using tent map, and extract symbol dynamics
	char nbits[50];
	for (int i=0; i<50; i++)
	{
		if (0.5 <= y)
		{
			y = 1.0 - y;
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

// the logistic compressor
double gist_cpr(double y, double K)
{
	// Iterate on y using the logistic map, and extract symbol dynamics
	char nbits[50];
	for (int i=0; i<50; i++)
	{
		if (0.5 <= y)
		{
			nbits[i] = 1;
		}
		else nbits[i] = 0;
		y *= 4.0*K*(1.0-y);
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

double xeig(double y, double K)
{
	if (1.0 < y) return 0.0;
	if (K < y) return 0.0;
	if (K-0.5 < y && y < 0.5) {
		return xeig(y*2.0*K, K);
	}
	return cpr(y,K) - 0.5;
}

double eig(double y, double K)
{
// y *= 3.0;
	y -= floor(y);
// if (K < y) return 0.0;
// if (0.5*K < y && y < 0.5) return 0.0;
// if (K < y) y -= 1;
// if (K < y) y -= K;
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
	double rng = atof(argv[2]);

#if 1
	int npts = 3803;
	for (int i=0; i<npts; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) npts);
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
		double ska = x / (2.0*lam);
		double skalo = cpr (ska, lam);
		double skahi = cpr (0.5+ska, lam);
		double lo = cpr (0.5*x, lam);
		double hi = cpr (0.5+0.5*x, lam);
		double skam = cpr (lam*x, lam);
		printf("%d	%g	%g	%g	%g	%g	%g	%g\n", i, x, y, skalo, skahi, lo, hi, skam);
#endif
#define CPR_SIM_EXT 1
#ifdef CPR_SIM_EXT
		x *= rng;
		double xx = x;
		// double y = cpr (x, lam);
		double y = extcpr (x, lam);
		double ey = extcpr (xx, lam);
		double sam = extcpr (xx/(2.0*lam), lam);
		printf("%d	%g	%g	%g	%g\n", i, x, y, ey, sam);
#endif
// #define EIG
#ifdef EIG
		double y = eig (x, lam);
		double lo = eig (x/(2.0*lam), lam);
		double hi = eig ((lam+x)/(2.0*lam), lam);
		printf("%d	%g	%g	%g	%g\n", i, x, y, lo, hi);
#endif

#ifdef ADJOINT
		// double y = tent_cpr (x, lam);
		// double z = gist_cpr (x, lam);
		double y = der (x, lam);
		double t = cpr(x, lam);
		double z = pdr(t, lam);
		printf("%d	%g	%g	%g\n", i, x, y, z);
#endif
	}
#endif

#ifdef LENGTH_OF_ALL_ONES
	// Verify the length formula for continguous runs of all-ones.
	int npts = 180;
	for (int i=0; i<npts; i++)
	{
		double lam = 0.5 + 0.5* (((double) i) + 0.5)/ ((double) npts);
		double w = len(1.0, lam);
		double p = flen(lam);
		p = ceil(p);
		printf("%d	%g	%g	%g %g\n", i, lam, w, p, w-p);
	}
#endif
}
