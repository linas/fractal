
/*
 * Distribution of gpf.
 * Linas June 2016
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gpf.h>

#define NBINS 20000
int cnt[NBINS];

void accum(int max)
{
	for (int i=0; i<NBINS; i++) cnt[i] = 0;

	for (int n=1; n<max; n++)
	{
		int g = gpf(n);
		double x = ((double) g) / ((double) n);
		x *= NBINS-1;
		int j = floor(x);
		cnt[j] ++;
	}
}

int
main(int argc, char* argv[])
{
	int m = atoi(argv[1]);
	accum(m);

	double norm = cnt[NBINS-1];
	double scale = 1.0 / NBINS;
	for (int i=1; i<NBINS; i++)
	{
		// printf("%d\t%d\n", i, cnt[i]);
		if (0 != cnt[i]) {
			double v = cnt[i] / norm;
			printf("%d\t%g\t0\tinf\n", i-1, (i-1)*scale);
			printf("%d\t%g\t%g\t%g\n", i, i*scale, v, 1.0/v);
			printf("%d\t%g\t0\tinf\n", i+1, (i+1)*scale);
		}
	}
}
