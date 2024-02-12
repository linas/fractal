/*
 * reforce.c
 * Build coherent states.
 * Unlike later work, these are not the piecewise coherent states,
 * and so they don't actually work out. See repiece.c for the
 * piece-wise variant.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double line (double x)
{
	return x-0.5;
}

// Recursively defined coherent function.
// Fast, easy, straight.
double coh(double beta, double omega, double y, int depth)
{
	if (0.5 * beta < y) return 0.0;
	if (0 == depth)
		return line(y/(0.5*beta));

	double xlo = y / beta;
	double xhi = xlo + 0.5;

	double dlo = coh(beta, omega, xlo, depth-1);
	double dhi = coh(beta, omega, xhi, depth-1);

	double ellie = (dlo + dhi) * omega;
	ellie += line(y/(0.5*beta));
	return ellie;
}

#ifdef USELESS_EQUIVALENT
// This function is similar to above, but just apply ell to
// recursive depth. Combined with coher sum below, this gives
// exactly the same result as coh() above. It's just a lot
// slower.
double ellen(double beta, double omega, double y, int depth)
{
	if (0.5 * beta < y) return 0.0;
	if (0 == depth)
		return line(y/(0.5*beta));

	double xlo = y / beta;
	double xhi = xlo + 0.5;

	double dlo = ellen(beta, omega, xlo, depth-1);
	double dhi = ellen(beta, omega, xhi, depth-1);

	double ellie = (dlo + dhi) * omega;
	return ellie;
}

// Perform explicit sum adding up recursive elts.
double coher(double beta, double omega, double y)
{
	double omen = 1.0;
	double sum = 0.0;
	for (int i=0; i<1000; i++)
	{
		sum += ellen(beta, omega, y, i);
		omen *= omega;
		if (omen < 1e-6) break;
	}
	return sum;
}
#endif

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s beta depth omega\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	int depth = atoi(argv[2]);
	double omega = atof(argv[3]);

	printf("#\n# beta=%g\n#\n", beta);

	omega /= beta;

#define NPTS 2019
	for (int j=0; j< NPTS; j++)
	{
		double x = (((double) j) + 0.5) / ((double) NPTS);
		double y = coh(beta, omega, x, depth);
		// double z = coher(beta, omega, x);
		printf("%d	%g	%g\n", j, x, y);
		fflush(stdout);
	}
}
