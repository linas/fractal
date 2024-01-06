/*
 * orbit.c
 *
 * Manual exploration of individual orbits.
 * Jan 2024
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#include "../periodic/selfie.c"

double teen(double beta, int k)
{
	double mid = 0.5* beta;
	for (int i=0; i<k; i++)
	{
		if (0.5 < mid) mid -= 0.5;
		mid *= beta;
	}
	return mid;
}

int cmp(const void* ida, const void* idb, void* pb)
{
	double beta = *((double*) pb);
	double ma = teen(beta, *((int*) ida));
	double mb = teen(beta, *((int*) idb));
	return (ma < mb) ? -1: 1;
}

void orbit(double beta, int order)
{
	double mid = 0.5* beta;
	for (int i=0; i<order; i++)
	{
		printf("%d	%g\n", i, mid);
		if (0.5 < mid) mid -= 0.5;
		mid *= beta;
	}
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
	orbit(beta, order);

	int mpts[20];
	for (int i=0; i< order; i++) mpts[i] = i;

	qsort_r(mpts, order, sizeof(int), cmp, &beta);

	printf("Midpoints:");
	for (int i=0; i<order; i++) printf(" %d", mpts[i]);
	printf ("\n");
}
