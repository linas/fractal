/*
 * ext-misc.C
 * Misc stuff.
 *
 * Linas Vepstas Oct 2020
 */

#include <math.h>
#include <stdio.h>

// See sidorov-bug.C for explanation of emrun
int emrun(double K)
{
   double beta = 2.0*K;
   double gold = 0.5 * (1.0 + sqrt(5));
   if (beta <= gold) return 1;

   double loga = (beta - 1.0) / (2.0-beta);
   loga = log(loga) / log(beta);
   loga = floor(loga) + 1.0;
   return (int) loga;
}

double tee(double x, double beta)
{
	if (x<0.5) return beta*x;
	return beta*(x-0.5);
}

double tau(double x, double beta)
{
	int em = emrun(0.5*beta);
	double a = 0.5/beta;
	double b = a * (1.0 + pow(beta, -em));
	if (a < x and x < b) return beta*x;
	if (a+0.5 < x and x < b+0.5) return beta*x;
	if (a+1.0 < x and x < b+1.0) return beta*x;
	if (x<0.5) return beta*x;

	if (x < 1.0) return beta*(x-0.5);
	return beta*(x-1.0);
}


int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K\n", argv[0]);
		exit (1);
	}
	double Kay = atof(argv[1]);

	int em = emrun(Kay);
	printf("#\n# K=%g m=%d\n#\n", Kay, em);

#define NBINS 501
#define SCALE 1.3333
	for (int i=0; i<NBINS; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) NBINS);
		x *= SCALE;
		double y = tau(x, 2.0*Kay);
		printf("%d	%g	%g\n", i, x, y);
	}
}

// ================================================================
