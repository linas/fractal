/*
 * irred-tree.c
 *
 * Create a graph of the bracket tree.
 *
 * December 2023
 */

#include "selfie.c"
#include "selfie-rational.c"
#include "selfie-tree.c"

// Given a position in a binary tree, return the corresponding beta
// index. This is the "good index" function, mapping the full binary
// tree to the trimmed tree of good indexes.
// Obsolete. Use good_index_map() in ndew code.
long moves_to_idx(long bitseq, int len)
{
	if (0 == bitseq%2)
	{
		printf("Error: expecting bitseq to be one-terminated (odd number)\n");
		return -2;
	}
	// printf("Enter frac = %ld / %d\n", bitseq, 1<<len);
	bitseq >>= 1;
	bitseq |= 1UL << len;
	return good_index_map(bitseq);
}

// Return the bit-sequence that approximates the rational p/q,
// truncated to len bits. Rationals have infinite periodic bitseqs,
// so this just expands that bitseq, and then truncates it.
//
// Bit strings are read left-to-right, with decimal point at far
// left. Bit strings interpreted as dyadics in the conventional sense.
// For example, 1/4 len=2 -> 0.01 and  3/8 len=3 -> 0.011
long rational_to_bracket_bitseq(int p, int q, int len)
{
	long bitseq = 0;
	for (int i=0; i<len; i++)
	{
		p <<= 1;
		if (q <= p) bitseq |= 1UL << (len-1-i);
		if (q <= p) p -= q;
	}
	return bitseq;
}

// Reverse the order of the bits in the bitseq
long bitseq_reverse(long bits, int len)
{
	long rev = 0;
	for (int i = 0; i<len; i++)
	{
		rev <<= 1;
		rev |= bits & 1UL;
		bits >>= 1;
	}
	return rev;
}

// Do the minkowski run-length encoding trick.
// Convert the bitseq to a continued fraction.
double bitseq_to_cf(long bitseq, int len)
{
	// printf("Enter cf frac = %ld / %d\n", bitseq, 1<<len);

	double frac = 0.0;

	int tail = bitseq & 0x1L;
	int bit = 0;
	// printf("bit %d\n", tail);
	int cnt = 1;
	for (int i=1; i<len; i++)
	{
		bit = (bitseq >> i) & 0x1L;
		// printf("bit %d\n", bit);
		if (bit == tail) cnt++;
		else
		{
			tail = bit;
			frac = 1.0 / ((double) cnt + frac);
			// printf("flip at bit %d frac = %g\n", i, frac);
			cnt = 1;
		}
	}
	frac = 1.0 / ((double) cnt+1 + frac);
	if (bit) frac = 1.0 - frac;
	// printf ("final frac = %g\n", frac);
	// printf("------\n");
	return frac;
}

int main(int argc, char* argv[])
{
// #define SPOT_CHECK
#ifdef SPOT_CHECK
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s <index>\n", argv[0]);
		exit(1);
	}
	long idx = atol(argv[1]);
	malloc_gold(idx+1);

	int len = 0;
	long bits = idx_to_moves(idx, &len);
	printf("Index=%ld ", idx);
	print_moves(bits, "bitseq=(", ")\n");

	long ridx = moves_to_idx(bits, len);
	printf("reconstructed index = %ld\n", ridx);
	if (idx != ridx)
		printf("Error: failed reconstruction!\n");

	double gold = find_gold(idx);
	printf("Index %ld  beta=%20.16g\n", idx, gold);
#endif

#define RATIONAL_EXPLORE
#ifdef RATIONAL_EXPLORE
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s <p> <q> <max-len>\n", argv[0]);
		exit(1);
	}
	int p = atoi(argv[1]);
	int q = atoi(argv[2]);
	int len = atoi(argv[3]);

	long maxidx = 1UL << (len+2);
	maxidx = 1UL << (len+6);
	//maxidx = 1UL << 34;
	//maxidx = 1UL << 30;
	maxidx = 1UL << 24;
	printf("Max idx = %ld\n", maxidx);

	printf("rational = %d/%d len=%d\n", p, q, len);

	long bits = rational_to_bracket_bitseq(p, q, len);

	// Terminal bit must always end in 1, so that bits is odd.
	bits |= 1UL;
	print_moves(bits, "Bracketing bits=(", ")\n");
	long idx = moves_to_idx(bits, len);
	double gold = golden_beta(idx);
	printf("Bracket index %ld  beta=%20.16g\n", idx, gold);

	printf("---\n");
	long orbits = rational_to_dyadic(p, q, len);
	orbits |= 1UL;
	print_moves(orbits, "orbit bits=(", ")\n");

	// Orbits MUST start with a leading one-bit.
	long lorbits = orbits;
	lorbits |= 1<<len;
	lorbits >>= 1;
	print_moves(lorbits, "leading-one bitseq=(", ")\n");

	long oridx = moves_to_idx(orbits, len);
	double orgold = golden_beta(oridx);
	printf("Orbit Index %ld  beta=%20.16g\n", oridx, orgold);
#endif

// #define VERIFY
#ifdef VERIFY
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s <max-order>\n", argv[0]);
		exit(1);
	}
	int maxord = atoi(argv[1]);
	long nmax = 1UL << maxord;
	malloc_gold(nmax);

	printf("Verify correctness to max order = %d\n", maxord);
	for (long idx=1; idx < nmax; idx++)
	{
		if (false == valid_gold_index(idx)) continue;

		int len = 0;
		long bits = idx_to_moves(idx, &len);
		// print_moves(bits, " # bits=(", ")\n");

		long ridx = moves_to_idx(bits, len);
		// printf("reconstruct %ld\n", ridx);

		if (ridx != idx)
		{
			printf("Error: bad indexing at idx=%ld rec=%ld", idx, ridx);
			print_moves(bits, " # bits=(", ")\n");
		}
	}
#endif

// #define PRINT_DYADIC_TO_BETA_MAP
#ifdef PRINT_DYADIC_TO_BETA_MAP
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s <max-order> <dyad-base>\n", argv[0]);
		exit(1);
	}
	int maxord = atoi(argv[1]);
	int dyad = atoi(argv[2]);

	long nmax = 1UL << maxord;

	printf("#\n# Bracket tree. max order = %d\n#\n", maxord);

	int overflo = 0;
	int maxdy = 1 << dyad;
	for (int i=1; i<maxdy; i++)
	{
		double x = ((double) i) / (double) maxdy;

		// Get the dyad
		int dlen = dyad;
		int bits = i;
		while (0 == bits %2)
		{
			bits >>= 1;
			dlen --;
		}
		long idx = moves_to_idx(bits, dlen);
		if (idx < 0) { overflo++; continue; }
		if (nmax <= idx) { overflo++; continue; }

		double gold = golden_beta(idx);

		// And now, a shifted version
		int shbit = bits;
		long loidx = moves_to_idx(shbit, dlen+1);
		if (loidx < 0) { overflo++; continue; }
		if (nmax <= loidx) { overflo++; continue; }
		double logold = find_gold(loidx);

		shbit = bits | 1<<dlen;
		long hiidx = moves_to_idx(shbit, dlen+1);
		double higold = find_gold(hiidx);

		double cf = bitseq_to_cf(bits, dlen);
		printf("%d	%d	%d	%ld	%g	%g	%g", i, bits, 1<<dlen, idx, x, gold, cf);
		printf("	%g", logold);
		printf("	%g", higold);

		print_moves(bits, " # bits=(", ")\n");
	}

	printf("#\n# Num overflows = %d\n#\n", overflo);
#endif
}

/* --------------------------- END OF LIFE ------------------------- */
