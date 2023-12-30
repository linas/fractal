/*
 * measure.c
 *
 * Explore sums over combs. Again.
 * A lot like bracket-finite.c but with dyadic fractions, instead.
 *
 * December 2023
 */

#include "selfie.c"
#include "selfie-util.c"
#include "selfie-tree.c"

// Sum over theta, starting from one.
// Very inefficient, but so what.
unsigned long comb_sum(unsigned long idx)
{
	unsigned long sum = 0UL;
	for (unsigned long n=1; n <= idx; n++)
	{
		bool ok = valid_gold_index(n);
		sum += ok;
	}
	return sum;
}

int main(int argc, char* argv[])
{
	// Map dyadic fractions to sums.
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <order>\n", argv[0]);
		exit(1);
	}
	int order = atoi(argv[1]);

	printf("#\n# Front dyadics. Order = %d\n#\n", order);
	int maxdy = 1 << order;
	for (int i=1; i<maxdy; i++)
	{
		double x = ((double) i) / (double) maxdy;

		int p = i;
		int q = 1UL << order;
		int gcf = gcd(p,q);
		p /= gcf;
		q /= gcf;

		// Convert dyadic to canonical tree numbering.
		int mcanon = (p + q) >> 1;

		unsigned long csum = comb_sum(mcanon);

		//unsigned long idx = good_index_map(mcanon);
		//if (MAXIDX < idx) continue;

		//unsigned long numo=0, deno=0;
		//tree_idx_to_dyafrac(idx, &numo, &deno);

		printf("%d	%g	%ld\n", i, x, csum);
	}
}
