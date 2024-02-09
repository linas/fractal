/*
 * bracket-finite.c
 *
 * Explore relationship between bracketing and mid-point moves.
 *
 * December 2023
 */

#include "selfie.c"
#include "selfie-util.c"
#include "selfie-rational.c"
#include "selfie-tree.c"

#include "necklace.h"

int main(int argc, char* argv[])
{
// #define MANUAL_EXPLORER
#ifdef MANUAL_EXPLORER
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <move>\n", argv[0]);
		exit(1);
	}
	long moves = atol(argv[1]);

	printf("Moves from index zero:\n");
	long down = good_index_map(moves);
	print_gold_info(down);

	printf("Just the index:\n");
	print_gold_info(moves);
#endif

// #define PRINT_PSI
#ifdef PRINT_PSI

	if (2 != argc) {
		fprintf(stderr, "Usage: %s <maxord>\n", argv[0]);
		exit(1);
	}
	int maxord = atoi(argv[1]);
	int prev = 0;
	for (int k=1; k<maxord; k++)
	{
		unsigned long mstart = 1UL << (k-1);
		unsigned long mend = 1UL << k;
		int rank = k+1;
		printf("Rank %d: ", rank);
		int tot = 0;
		for (unsigned long m=mstart; m<mend; m++)
		{
			if (false == valid_gold_index(m)) continue;
			printf("%ld, ", m);
			// printf("%ld, ", m-prev);
			// if (1 < m-prev) printf("%ld, ", m-prev);
			prev = m;
			tot++;
		}
		printf("\n");
		int moreau_nu = necklace(rank);
		printf("Above has %d entries; expected %d\n\n", tot, moreau_nu);
		if (0 != moreau_nu - tot) printf("Error!!!!!!! XXXXXXXXXXXX\n");
	}
#endif

#define PRINT_BFILE
#ifdef PRINT_BFILE
	// To generate OEIS Bfile.
	printf("# Generated with https://github.com/linas/fractal/blob/master/experiments/beta-xform/periodic/bracket-finite.c\n");
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <nterms>\n", argv[0]);
		exit(1);
	}
	int nterms = atoi(argv[1]);
	int tot = 0;
	long idx = -1;
	while (tot <= nterms)
	{
		idx ++;
		if (false == valid_gold_index(idx)) continue;
		printf("%d %ld\n", tot, idx);
		tot++;
	}
	printf("\n");
#endif

// #define INDICATOR_GRAPH
#ifdef INDICATOR_GRAPH
	// Print the indicator function, as well as the front values.
	// This is used for graphs in the paper.
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <maxord>\n", argv[0]);
		exit(1);
	}
	int maxord = atoi(argv[1]);
	double fact = 1.0;
	for (int k=0; k<maxord; k++)
	// for (int k=maxord-1; k<maxord; k++)
	{
		fact *= k+1;
		long comb_sum = 0;
		double idxsum = 0.0;

		unsigned long mstart = 1UL << k;
		unsigned long mend = 1UL << (k+1);

		// First loop, for normalization
		for (unsigned long m=mstart; m<mend; m++)
		{
			bool ok = valid_gold_index(m);
			comb_sum += ok;
			long idx = good_index_map(m);
			idxsum += idx;
		}
		// printf("# nu=%d mstart=%ld comb=%ld idxs=%ld\n", k+2, mstart, comb_sum, idxsum);
		printf("%d	%ld	%g	%ld	%g\n", k+2, mstart, fact, comb_sum, idxsum);

#if 0
		// int nu = k+2
		// Comb_norm is Moreau necklace at rank nu = (k+2)
		// So k=0 |-> nu=2 |-> Moreau = 1
		double comb_norm = comb_sum;
		comb_sum = 0;

		// Factor of four because we want 2^nu and mstart = 2^(nu-2)
		// Note that expo goes to  1/nu for nu larger than about 15
		double expo = comb_norm / (4.0*mstart);
		// double expo = 1.0 / ((double) (k+2));

		double inorm = idxsum;
		idxsum = 0.0;

		for (unsigned long m=mstart; m<mend; m++)
		{
			bool ok = valid_gold_index(m);
			comb_sum += ok;
			long idx = good_index_map(m);
			idxsum += idx;

			double comb_meas = ((double) comb_sum) / comb_norm;
			double imeas = idxsum / inorm;

			double ska = pow(comb_meas, expo);
			double ex = (((double) (m-mstart))+0.5) / ((double)(mend-mstart));

			printf("%ld	%d	%ld	%ld	%g	%g	%g	%g	%g\n",
				m, ok, comb_sum, idx, idxsum, comb_meas, ska, ex, imeas);
		}
		printf("\n");
#endif
		fflush(stdout);
	}
#endif

// #define FIND_BIGGEST_FRONT
#ifdef FIND_BIGGEST_FRONT
	// Find the largest value for the front in a given rank.
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <maxord>\n", argv[0]);
		exit(1);
	}
	int maxord = atoi(argv[1]);

	for (int k=1; k<maxord; k++)
	{
		unsigned long maxlead = 0;
		unsigned long mov = 0;

		unsigned long start = 1UL << k;
		unsigned long end = 1UL << (k+1);
		for (unsigned long m=start; m<end; m++)
		{
			unsigned long idx = good_index_map(m);
			if (maxlead < idx)
			{
				mov = m;
				maxlead = idx;
			}
		}

		// printf("order = %d mov= %ld max= %ld\n", k, mov, maxlead);
		printf("%d	%ld	%ld\n", k, mov, maxlead);
	}
#endif

#ifdef BETA_RANKS
	// Dump the gold values for front values into a file.
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <maxord>\n", argv[0]);
		exit(1);
	}
	int maxord = atoi(argv[1]);

	for (int k=1; k<maxord; k++)
	{
		unsigned long start = 1UL << k;
		unsigned long end = 1UL << (k+1);
		for (unsigned long m=start; m<end; m++)
		{
			unsigned long idx = good_index_map(m);
			if (MAXIDX < idx) continue;
			double gold = golden_beta(idx);

			printf("%ld	%d	%ld	%g\n", m, k, idx, gold);
		}
	}
#endif

// #define UNIFORM_CONVERGENCE
#ifdef UNIFORM_CONVERGENCE
	// Look at how the gold sequence converges.
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <maxord>\n", argv[0]);
		exit(1);
	}
	int maxord = atoi(argv[1]);

	unsigned long start = 1UL << maxord;
	unsigned long end = 1UL << (maxord+1);
	double delta = end-start;
	delta = 1.0 / delta;
	for (unsigned long m=start; m<end; m++)
	{
		double x = (((double) m-start) +0.5) * delta;
		unsigned long idx = good_index_map(m);
		if (MAXIDX < idx) continue;
		double gold = golden_beta(idx);
		unsigned long lidx = bracket_gold_left(idx);
		unsigned long ridx = bracket_gold_right(idx);
		double lgold = golden_beta(lidx);
		double rgold = golden_beta(ridx);
		double diff = rgold - lgold;
		printf("%ld	%g	%g	%g\n", m, x, gold, diff);
#ifdef TOO_COMPLICATED
		printf("%ld	%g	%16.14g", m, x, gold);
		for (int k=1; k<8; k++)
		{
			// idx = bracket_gold_left(idx);
			idx = bracket_gold_right(idx);
			gold = golden_beta(idx);
			printf("	%16.14g", gold);
		}
		printf("\n");
#endif
	}
#endif

// #define WORST_CONVERGENT
#ifdef WORST_CONVERGENT
	// Explore what is happening at beta=1
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <maxord>\n", argv[0]);
		exit(1);
	}
	int maxord = atoi(argv[1]);

	for (int k=1; k<maxord; k++)
	{
		unsigned long left = 1UL << k;
		unsigned long idx = good_index_map(left);
		double gold = golden_beta(idx);
		unsigned long ridx = good_index_map(left+1);
		double rgold = golden_beta(ridx);

		printf("%d	%16.12g	%16.12g	%g\n", k, gold, rgold, rgold-gold);
	}
#endif
}
