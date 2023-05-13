
/*
 * Sampling in the transfer operator.
 *
 * Dec 2017 Linas Vepstas
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NBINS 979
double cnt[NBINS];

double xiter (double x, double K, int lvl)
{
	if (K < x)
	{
		int bin = floor(x*NBINS);
		cnt[bin] += 1.0;
		return 0.0;
	}
	if (lvl < 0)
	{
		int bin = floor(x*NBINS);
		cnt[bin] += 1.0;

		return 1.0/K;
	}
	double otk = 0.5 / K;
	double xtk = otk*x;
	lvl --;
	double sum = 0.0;
	sum += xiter(xtk, K, lvl);
	sum += xiter(xtk+0.5, K, lvl);
	sum *= otk;
	return sum;
}

int main(int argc, char* argv[])
{
	double K = atof(argv[1]);
	double xmin = atof(argv[2]);
	double xmax = atof(argv[3]);
	double nsamp = atof(argv[4]);

	for (int i=0; i< NBINS; i++) cnt[i] = 0.0;

#define NREC 26
	for (int i=0; i<nsamp; i++)
	{
		double x = rand();
		x /= RAND_MAX;
		x *= (xmax - xmin);
		x += xmin;

		xiter(x, K, NREC);
	}
	for (int i=0; i< NBINS; i++)
	{
		cnt[i] *= NBINS / nsamp;
		double x = ((double) i + 0.5) / ((double) NBINS);
		printf ("%d	%g	%g", i, x, cnt[i]);
		printf ("\n");
	}
}
