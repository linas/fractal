/*
 * irred-allowed.c
 *
 * Verify it allowed-sequence bracketing formulas.
 *
 * December 2023
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "irred-gold.c"

// =================================================================

// Compute the period (order) of the iteration of the midpoint
// This is just sanity check; beta should "already be correct".
int midpoint_period(double beta)
{
#define DELTA 1e-15
#define EPSI  1e-13
	int count = 0;
	double midpoint = 0.5-DELTA;
	do
	{
		count ++;
		midpoint = tee(beta, midpoint);
		// printf("%d  %g\n", count, midpoint);
	}
	while (fabs(midpoint-0.5)>EPSI && count < 60);
	return count;
}

// =================================================================

// Compute the bitsequence of the iteration of the midpoint
long midpoint_bitseq(double beta)
{
#define DELTA 1e-15
#define EPSI  1e-13
	long bitseq = 0;
	int count = 0;
	double midpoint = 0.5-DELTA;
	do
	{
		midpoint = tee(beta, midpoint);
		// printf("%d  %g\n", count, midpoint);
		if (midpoint < 0.5) bitseq |= 1UL;
		bitseq <<= 1;
		count ++;
	}
	while (fabs(midpoint-0.5)>EPSI && count < 60);
	return bitseq;
}

// =================================================================

int main(int argc, char* argv[])
{
	// Obtain index from command line.
	if (1 == argc) {
		fprintf(stderr, "Usage: %s <index>\n", argv[0]);
		exit(1);
	}

	long idx = atol(argv[1]);

	malloc_gold(100);
	double gold = find_gold(idx);
	printf("Index: %ld gold=%g", idx, gold);

	int orbitp = midpoint_period(gold);
	int ord = order(idx);
	printf("midpoing period=%d order=%d\n", orbitp, ord);

	long bitseq = midpoint_bitseq(gold);

	prt_bitstr(bitseq, "bitseq ", "\n");
}
