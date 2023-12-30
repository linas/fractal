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

#include "necklace.h"

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

// Reset sum to one at every order.
// This is S_\nu from the paper, wiht \nu infered from the idx.
unsigned long comb_sum_order(unsigned long idx)
{
	int order = bitlen(idx);
	unsigned long nstart = 1UL << (order-1);

	unsigned long sum = 0UL;
	for (unsigned long n=nstart; n <= idx; n++)
	{
		bool ok = valid_gold_index(n);
		sum += ok;
	}
	return sum;
}

// A_\nu from the paper
double comb_sum_norm(unsigned long idx)
{
	unsigned long sum = comb_sum_order(idx);

	int rank = bitlen(idx) + 1;
	double moreau_nu = (double) necklace(rank);

	double renorm = ((double) sum) / moreau_nu;

	double scale = moreau_nu / ((double) (1UL << rank));
	double rescale = pow (renorm, scale);

	return rescale;
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

		// unsigned long csum = comb_sum(mcanon);
		// unsigned long csum = comb_sum_order(mcanon);
		double nsum = comb_sum_norm(mcanon);

		//unsigned long idx = good_index_map(mcanon);
		//if (MAXIDX < idx) continue;

#if 0
		unsigned long numo=0, deno=0;
		tree_idx_to_dyafrac(csum, &numo, &deno);
		double csum_dya = ((double) numo) / ((double) deno);
#endif

		printf("%d	%g	%g\n", i, x, nsum);
	}
}
