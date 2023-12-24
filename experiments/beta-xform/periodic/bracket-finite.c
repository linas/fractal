/*
 * bracket-finite.c
 *
 * Explore relationship between bracketing and mid-point moves.
 *
 * December 2023
 */

#include "selfie.c"

int main(int argc, char* argv[])
{
#ifdef MANUAL_EXPLORER
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <move>\n", argv[0]);
		exit(1);
	}
	long moves = atol(argv[1]);

	printf("Moves from index zero:\n");
	long down = move_gold_index(moves);
	print_gold_info(down);

	printf("Just the index:\n");
	print_gold_info(moves);
#endif

// #define GRAPH_FOR_PAPER
#ifdef GRAPH_FOR_PAPER
	// Dump the front values into a file, for graphing. Used for
	// graphs in the paper.
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <maxord>\n", argv[0]);
		exit(1);
	}
	int maxord = atoi(argv[1]);
	long maxind = 1UL << maxord;

	long idxsum = 0;
	long tsum = 0;
	int ord = 1;
	for (long mov=1; mov< maxind; mov++)
	{
		if (mov%ord == 0)
		{
			ord <<= 1;
			tsum = 0;
		}
		bool ok = valid_gold_index(mov);
		tsum += ok;
		long idx = front_sequence(mov);
		idxsum += idx;

		printf("%ld	%d	%ld	%ld	%ld\n", mov, ok, tsum, idx, idxsum);
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
			unsigned long idx = front_sequence(m);
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
			unsigned long idx = front_sequence(m);
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
		unsigned long idx = front_sequence(m);
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

#define WORST_CONVERGENT
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
		unsigned long idx = front_sequence(left);
		double gold = golden_beta(idx);
		unsigned long ridx = front_sequence(left+1);
		double rgold = golden_beta(ridx);

		printf("%d	%16.12g	%16.12g	%g\n", k, gold, rgold, rgold-gold);
	}
#endif
}
