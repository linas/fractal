/*
 * irred-midpoint.c
 *
 * Verify the orbit exapansions and orbit labeling formulas.
 * The mid-point should orbit to the left and right exactly the same
 * way as the bitsequence that defines the polynomial for that beta.
 * Make sure everything works as claimed.
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
	long bitseq = 1;
	int count = 0;
	double midpoint = 0.5*beta-DELTA;
	do
	{
		bitseq <<= 1;
		midpoint = tee(beta, midpoint);
		if (midpoint > 0.5-EPSI) bitseq |= 1UL;
		// printf("%d  %7.5g lo=%d 0x%lx\n", count, midpoint, (midpoint > 0.5-EPSI), bitseq);
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
	bool valid = is_valid_index(idx);
	double pzero = find_poly_zero(idx);
	printf("Index: %ld valid=%d pzero=%g\n", idx, valid, pzero);

	if (false == valid)
	{
		long factor = theta_factor(idx, pzero);
		printf("Invalid index, factors to %ld\n", factor);
	}

	int orbitp = midpoint_period(pzero);
	int ord = order(idx);
	printf("midpoint period=%d order=%d\n", orbitp, ord);

	long bitseq = midpoint_bitseq(pzero);
	int bord = order(bitseq/2 -1);
	printf("Orbit bitseq = 0x%lx order=%d ", bitseq, bord);
	for (int i=0; i<bord; i++)
		printf("%ld", 0x1 & bitseq>>(bord-i-1));
	printf("\n");

	printf("Orbit is ");
	print_orbit(bord+1, pzero);

	printf("yes, its %ld\n", bitseq/2);
}
