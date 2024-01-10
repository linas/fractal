/*
 * visitation.c
 *
 * Visitation function, from the diary.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "visits.c"

#define DEPTH 5000
unsigned long visits[DEPTH];
double midpts[DEPTH];
int sorted[DEPTH];

static int cmp(const void* ida, const void* idb)
{
	int ia = *((int*) ida);
	int ib = *((int*) idb);
	return (midpts[ia] < midpts[ib]) ? -1: 1;
}

/*
 * Arrange midpoint iterations into a tree, so that each midpoint
 * splits a prior midpoint interval into two. This generates a kind
 * of prefered Borel set, where the boundaries are always at midpoints.
 */
void vtree(double beta)
{
	int tot = visitation_tree(beta, visits, midpts, DEPTH, -1.0);

	// Sort into sequential order.
	for (int i=1; i<DEPTH; i++) sorted[i] = i;

	qsort(sorted, tot, sizeof(int), cmp);
}

#ifndef NOMAIN
int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s beta\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);

	vtree(beta);

	printf("%d	%d	%ld	%g	%g\n", -1, -1, 0UL, 0.0, 0.0);
	for (int i=0; i<DEPTH; i++)
	{
		if (NEG_ONE == visits[i]) break;
		int j = sorted[i];
		unsigned long v = visits[j];
		double dya = canonical_dyadic(v);
		double midpnt = midpts[j];
		printf("%d	%d	%ld	%g	%g\n", i, j, v, dya, midpnt);
	}
	printf("%d	%d	%ld	%g	%g\n", -1, -1, 0UL, 1.0, 1.0);

#ifdef STEPPY
	printf("%d	%d	%ld	%g	%g\n", -1, -1, 0UL, 0.0, 0.0);
	double yprev = 0.0;
	for (int i=0; i<DEPTH; i++)
	{
		if (NEG_ONE == visits[i]) break;
		int j = sorted[i];
		unsigned long v = visits[j];
		double dya = canonical_dyadic(v);
		printf("%d	%d	%ld	%g	%g\n", i, j, v, dya, yprev);
		printf("%d	%d	%ld	%g	%g\n", i, j, v, dya, midp[j]);
		yprev = midpts[j];
	}
	printf("%d	%d	%ld	%g	%g\n", -1, -1, 0UL, 1.0, 1.0);
#endif
}
#endif
