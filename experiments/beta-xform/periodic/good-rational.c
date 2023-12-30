/*
 * good-rational.c
 *
 * Explore the good-rational mapping, which maps rationals to those
 * with valid expansions.
 *
 * December 2023
 */

#include "selfie.c"
#include "selfie-rational.c"
#include "selfie-tree.c"

// Try to guess a rational.
void low_guess(double rat, int* p, int* q, int ndigits)
{
	int maxdeno = 1UL << ndigits;
	double epsi = 1.0 / ((double) maxdeno);
	for (int i=2; i<maxdeno; i++)
	{
		if ((fabs(fmod(rat*i, 1.0)) < epsi) ||
		    (fabs(1.0-fmod(rat*i, 1.0)) < epsi))
		{
			*q = i;
			*p = floor(rat*i + 1.5*epsi);
			return;
		}
	}
}

// This implements the good-rational map, that accepts any rational as
// input, and maps it to theta-bar, the set of acceptable rationals that
// have valid self-describing periodic orbits.
//
// The implementation is via approximation: the input rational is
// converted to a bit-string. This bitstring is interpreted as a
// sequence of moves down the binary tree. The moves are taken, using
// the leader function for the right-side moves. This is equivalent
// to passing the input bitstring as a canonical integer to the
// good-index map. If the good-index map overflows, then *quot=0
// and this is a conversion failure.
//
// Otherwise, good_index_map() returns a very large integer that can
// be interpreted as a sequence of moves down the trimmed tree. Taking
// those moves gives us some dyadic rational, with some huge denominator,
// that approximates the desired output. We now have to reconstruct
// what the intended limit of the approximant was. I suppose that a
// principled approach could be used, to take this limit. But we do
// a quick cheap hack: just look for the nearest rational that has some
// relatively small denominator. Yuck, but it works well-enough for
// graphs.

void good_rational_map(int pin, int qin, int* pout, int* qout)
{
	// Reduce to lowest terms
	unsigned long gcf = gcd(pin, qin);
	pin = pin / gcf;
	qin = qin / gcf;

	// Attempt to find a high-quality approximation.
	int minlen = 16;
	unsigned long idx = NEG_ONE;
	for (int len=40; minlen < len; len -=3)
	{
		unsigned long moves = rational_to_moves(pin, qin, len);
		idx = good_index_map(moves);

		// Once more, for the win.
		if (idx < MAXIDX)
		{
			moves = rational_to_moves(pin, qin, len-3);
			idx = good_index_map(moves);
			break;
		}
	}

	// Approximation failure.
	if (MAXIDX < idx)
	{
		*qout = 0;
		return;
	}

	// Convert the approximant back to a rational
	unsigned long tno = 2UL * idx + 1UL;
	unsigned long pp, qq;
	moves_to_rational(tno, &pp, &qq);
	double orat = ((double) pp) / ((double) qq);

	// Reduce the rational to lowest terms, given that
	// the approximant was probably off by a few bits.
	// It might be better to actually twiddle the approximant
	// directly, but this will do, for now.
	low_guess(orat, pout, qout, minlen-2);
}

int main(int argc, char* argv[])
{
// #define RATIONAL_EXPLORE
#ifdef RATIONAL_EXPLORE
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s <p> <q> <max-len>\n", argv[0]);
		exit(1);
	}
	int p = atoi(argv[1]);
	int q = atoi(argv[2]);
	int len = atoi(argv[3]);
	printf("\n");

	double rat = ((double) p) / ((double) q);
	printf("Input rational = %d/%d len=%d as float=%18.16g\n", p, q, len, rat);

	// Truncate the dyadics.
	unsigned long gcf = gcd(p, q);
	p /= gcf;
	q /= gcf;
	if (is_power_of_two(q))
		len = bitlen(q) - 2;

	long moves = rational_to_moves(p, q, len);
	print_moves(moves, "Tree moves=(", ")\n");

	unsigned long pp, qq;
	moves_to_rational(moves, &pp, &qq);
	double rrat = ((double) pp) / ((double) qq);
	printf("Reconstructed rational = %ld/%ld = %18.16g\n", pp, qq, rrat);

	unsigned long idx = good_index_map(moves);
	double gold = golden_beta(idx);
	printf("Bracket index %ld  beta=%20.16g\n", idx, gold);

	unsigned long oradic = beta_to_dyadic(gold);
	print_dyadic(oradic, 63, "Actual orbit: ", "\n");

	unsigned long oridx = beta_to_index(gold);
	printf("Reconstructed index from orbit: %ld\n", oridx);
	if (oridx != idx)
	{
		double rgold = golden_beta(oridx);
		printf("Warning: indexes don't match! diff=%g\n", gold-rgold);
	}

	unsigned long tno = 2UL * idx + 1UL;

	// The truncated expansion for dyadics hurts. So don't do that.
	if (is_power_of_two(q)) tno = idx;

	moves_to_rational(tno, &pp, &qq);
	double orat = ((double) pp) / ((double) qq);
	printf("Orbit rational = %ld/%ld = %18.16g\n", pp, qq, orat);

	int po=0, qo=0;
	low_guess(orat, &po, &qo, 14);
	printf("guess its orbit is %d/%d\n", po,qo);
	printf("\n");

#endif

#define GOOD_GRAPH
#ifdef GOOD_GRAPH
	// Implement above explorer as a direct map
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s <npts>\n", argv[0]);
		exit(1);
	}
	int npts = atoi(argv[1]);

	for (int i=0; i<= npts; i++)
	{
		// Input rational
		int pi = i;
		int qi = npts;

		// Output rational
		int po=0, qo=0;
		good_rational_map(pi, qi, &po, &qo);

		if (0 == qo) continue;

		double irat = ((double) pi) / ((double) qi);
		double orat = ((double) po) / ((double) qo);
		printf("%d	%d	%g	%d	%d	%g\n", pi, qi, irat, po, qo, orat);
		fflush (stdout);
	}
#endif
}

/* --------------------------- END OF LIFE ------------------------- */
