/*
 * orbit.c
 *
 * Manual exploration of individual orbits.
 * Jan 2024
 */

#include <stdio.h>
#include <stdlib.h>

#include "../periodic/selfie.c"

void orbit(double* mpts, double beta, int order)
{
	double mid = 0.5* beta;
	for (int i=0; i<order; i++)
	{
		printf("%d	%g\n", i, mid);
		mpts[i] = mid;
		if (0.5 < mid) mid -= 0.5;
		mid *= beta;
	}
}

int cmp(const void *ma, const void* mb)
{
	return *((double *) ma) < *((double *) mb);
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s idx\n", argv[0]);
		exit (1);
	}
	int idx = atoi(argv[1]);

	double beta = golden_beta(idx);
	printf("idx=%d beta=%g\n", idx, beta);

	int order = bitlen(idx) +1;
	double mpts[20];
	orbit(mpts, beta, order);

	qsort(mpts, sizeof(double), order, cmp);
}
