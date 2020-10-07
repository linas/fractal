/*
 * ext-iterate.C
 * Iterate the extended map
 *
 * Linas Vepstas Oct 2020
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "emrun.C"

double tee(double x, double beta)
{
	if (x<0.5) return beta*x;
	return beta*(x-0.5);
}

double tau(double x, double beta, double a, double b)
{
	if (a < x and x < b) return beta*x;
	if (a+0.5 < x and x < b+0.5) return beta*x;
	if (a+1.0 < x and x < b+1.0) return beta*x;
	if (x<0.5) return beta*x;

	if (x < 1.0) return beta*(x-0.5);
	return beta*(x-1.0);
}

#define NBINS 503
double histo[NBINS+1];

int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K\n", argv[0]);
		exit (1);
	}
	double Kay = atof(argv[1]);
	double beta = 2.0*Kay;

	int em = emrun(Kay);
	double a = 0.5/beta;
	double b = a * (1.0 + pow(beta, -em));
	printf("#\n# K=%g m=%d\n#\n", Kay, em);

	for (int i=0; i<=NBINS; i++)
		histo[i] = 0.0;

#define SCALE 1.333333
#define NSAMP 10000
	for (int i=0; i<NBINS*NSAMP; i++)
	{
		long int r = random();
		double x = ((double) r) / ((double) RAND_MAX);
		x *= SCALE;

		for (int j=0; j<50; j++)
		{
			x = tau(x, beta, a, b);

			double xb = NBINS * x / SCALE;
			int ibin = xb;
			histo[ibin] += 1.0;
		}
	}

	double acc = 0.0;
	for (int i=0; i<NBINS; i++)
		acc += histo[i];

	for (int i=0; i<NBINS; i++)
		histo[i] *= NBINS/acc;

	for (int i=0; i<NBINS; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) NBINS);
		x *= SCALE;
		double y = histo[i];
		printf("%d	%g	%g\n", i, x, y);
	}
}

// ================================================================
