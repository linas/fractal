/*
 * orbit.c
 *
 * Manual exploration of individual orbits.
 * Print out orbits, and then sort them into order.
 * Jan 2024
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>

#include "../periodic/selfie.c"

double teen(double beta, int k)
{
	if (k < 0) return 1.0;

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

	// Print the orbit
	int order = bitlen(idx) +1;
	orbit(beta, order);

	// Print midpoints in sorted order.
	int mpts[20];
	for (int i=0; i< order; i++) mpts[i] = i;
	qsort_r(mpts, order, sizeof(int), cmp, &beta);

	printf("Midpoints:");
	for (int i=0; i<order; i++) printf(" %d", mpts[i]);
	printf ("\n");

	double xlo = 0.0;
	for (int i=0; i< order; i++)
	{
		double xhi = teen(beta, mpts[i]);

		double leflo = xlo / beta;
		double lefhi = xhi / beta;
		double riglo = leflo + 0.5;
		double righi = lefhi + 0.5;
		printf("I %d [%g %g] \tL=[%g %g]  \tR=[%g %g]\n",
			order-i, xlo, xhi, leflo, lefhi, riglo, righi);

		xlo = xhi;
	}
}
