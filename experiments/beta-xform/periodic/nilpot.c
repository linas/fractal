
/*
 * nilpot.c
 *
 * Explore the niplotent irrationals.
 * January 2018.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s NMAX\n", argv[0]);
		exit (1);
	}
	int nmax = atoi(argv[1]);

#define NBINS 1501
	double cnt[NBINS];
	for (int i=0; i< NBINS; i++) cnt[i] = 0;

	for (int n=2; n< nmax; n++)
	{
		double on = 1.0 / ((double) n);
		int m = 0;
		while (1)
		{
			m++;
			double K = 0.5 * pow (2.0*m+1.0, on);
			// printf("its %d %d %g\n", n, m, K);
			if (1.0 <= K) break;
			int bin = K*NBINS;
			cnt[bin] += 1.0;
		}
	}
	double tot = 0;
	for (int i=0; i< NBINS; i++) tot += cnt[i];
	for (int i=0; i< NBINS; i++) cnt[i] /= tot;
	for (int i=0; i< NBINS; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NBINS);
		printf("%d	%g	%g\n", i, x, cnt[i]);
	}
}
