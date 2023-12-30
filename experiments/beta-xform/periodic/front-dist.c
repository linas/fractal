/*
 * front-dist.c
 *
 * Explore what the fronts are doing.
 * A lot like bracket-finite.c but with dyadic fractions, instead.
 * More like bracket-steps.c but with explicit dyadics.
 *
 * December 2023
 */

#include "selfie.c"
#include "selfie-util.c"
#include "selfie-tree.c"

int main(int argc, char* argv[])
{
	// Map front values to dyadics.
	// This is used for graphs in the paper.
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <order>\n", argv[0]);
		exit(1);
	}
	int order = atoi(argv[1]);

	double idxsum = 0.0;

	printf("#\n# Front dyadics. Order = %d\n#\n", order);
	int maxdy = 1 << order;
	// for (int i=1; i<maxdy; i++)
	for (int i=maxdy-1; 0 < i; i--)
	{
		double x = ((double) i) / (double) maxdy;

		int p = i;
		int q = 1UL << order;
		int gcf = gcd(p,q);
		p /= gcf;
		q /= gcf;

		// Convert dyadic to canonical tree numbering.
		int mcanon = (p + q) >> 1;

		unsigned long idx = good_index_map(mcanon);
		if (MAXIDX < idx) continue;

		idxsum += idx;
		double hei = log2(idx);
		double sei = log2(idxsum);

		//unsigned long numo=0, deno=0;
		//tree_idx_to_dyafrac(idx, &numo, &deno);

		printf("%d	%g	%ld	%g	%g\n", i, x, idx, hei, sei);
	}
}
