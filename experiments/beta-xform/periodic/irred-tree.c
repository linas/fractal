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

// cheap hack
void low_guess(double rat, int* p, int* q)
{
	for (int i=2; i<16386; i++)
	{
		if ((fabs(fmod(rat*i, 1.0)) < 1e-6) ||
		    (fabs(1.0-fmod(rat*i, 1.0)) < 1e-6))
		{
			*q = i;
			*p = floor(rat*i + 1e-5);
			return;
		}
	}
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

	long bits = idx_to_moves(idx);
	printf("Index=%ld ", idx);
	print_moves(bits, "Move sequence=(", ")\n");

	unsigned long p;
	unsigned long q;
	moves_to_rational(bits, &p, &q);
	printf("Matching rational = %ld / %ld\n", p, q);

	long ridx = good_index_map(bits);
	printf("reconstructed index = %ld\n", ridx);
	if (idx != ridx)
		printf("Error: failed reconstruction!\n");

	double gold = golden_beta(idx);
	printf("Index %ld  beta=%20.16g\n", idx, gold);
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

	printf("Verify correctness to max order = %d\n", maxord);
	for (long idx=1; idx < nmax; idx++)
	{
		if (false == valid_gold_index(idx)) continue;

		int len = 0;
		long bits = idx_to_moves(idx, &len);
		// print_moves(bits, " # bits=(", ")\n");

		bits >>= 1;
		bits |= 1UL << len;
		long ridx = good_index_map(bits);
		// printf("reconstruct %ld\n", ridx);

		if (ridx != idx)
		{
			printf("Error: bad indexing at idx=%ld rec=%ld", idx, ridx);
			print_moves(bits, " # bits=(", ")\n");
		}
	}
#endif

#define PRINT_DYADIC_TO_BETA_MAP
#ifdef PRINT_DYADIC_TO_BETA_MAP
	// Generate datafile for the bracket-tree graph, and the good-dyad map.
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s <order>\n", argv[0]);
		exit(1);
	}
	int order = atoi(argv[1]);

	printf("#\n# Bracket tree. Order = %d\n#\n", order);
	printf("# Columns:\n");
	printf("# i  p  q  canon-n good-idx p/q beta\n#\n");

	int overflo = 0;
	int maxdy = 1 << order;
	for (int i=1; i<maxdy; i++)
	{
		double x = ((double) i) / (double) maxdy;

		// Reduce fraction to lowest order. This will have
		// the form p / 2^n where p is odd.
		int n = order;
		int p = i;
		while (0 == p %2)
		{
			p >>= 1;
			n --;
		}

		// Convert dyadic to canonical tree numbering.
		// This is 1 at the root, 2,3 in first row, 4,5,6,7 in next row, etc.
		// This is m + 2^(n-1)  where p = 2m+1
		int mcanon = (p >> 1) + (1UL << (n-1));

		unsigned long idx = good_index_map(mcanon);
		if (MAXIDX < idx) { overflo++; continue; }

		double gold = golden_beta(idx);

		printf("%d	%d	%d	%d	%ld	%g	%g", i, p, 1<<n, mcanon, idx, x, gold);

		print_moves(mcanon, " # bits=(", ")\n");
	}

	printf("#\n# Num overflows = %d\n#\n", overflo);
#endif

// #define EXPERIMENTAL_DYADIC_TO_BETA_MAP
#ifdef EXPERIMENTAL_DYADIC_TO_BETA_MAP
	// Old code to generate datafile for the original bracket-tree graph.
	// Includes extra grunge for experiments. Don't use this; use the
	// new improved version above.
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s <order>\n", argv[0]);
		exit(1);
	}
	int order = atoi(argv[1]);

	printf("#\n# Bracket tree. Order = %d\n#\n", order);

	int overflo = 0;
	int maxdy = 1 << order;
	for (int i=1; i<maxdy; i++)
	{
		double x = ((double) i) / (double) maxdy;

		// Reduce fraction to lowest order.
		int dlen = order;
		int bits = i;
		while (0 == bits %2)
		{
			bits >>= 1;
			dlen --;
		}
		int shbit = bits;
		shbit >>= 1;
		shbit |= 1UL << dlen;
		unsigned long idx = good_index_map(shbit);
		if (MAXIDX < idx) { overflo++; continue; }

		double gold = golden_beta(idx);

		// And now, a shifted version
		shbit = bits;
		shbit >>= 1;
		shbit |= 1UL << (dlen+1);
		long loidx = good_index_map(bits);
		if (loidx < 0) { overflo++; continue; }
		double logold = golden_beta(loidx);

		shbit = bits | 1<<dlen;
		shbit >>= 1;
		shbit |= 1UL << (dlen+1);
		long hiidx = good_index_map(bits);
		double higold = golden_beta(hiidx);

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
