
/*
 * Distribution of gpf.
 * Linas June 2016
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gpf.h>

#define NBINS 2000
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

	for (int i=0; i<NBINS; i++)
	{
		printf("%d\t%d\n", i, cnt[i]);
	}
}
