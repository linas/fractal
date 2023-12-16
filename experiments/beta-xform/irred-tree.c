/*
 * irred-tree.c
 *
 * Create a graph of the bracket tree.
 *
 * December 2023
 */

#include "irred-gold.c"

// Give a sequence of 'len' left-right moves encoded as bits in bitseq
// return the corresponding beta index. Thus, bitseq is a dyadic rational
// of bitseq/2^len
long bitseq_to_idx(long bitseq, int len)
{
	long front = 1;
	// printf("Enter frac = %ld / %d\n", bitseq, 1<<len);
	for (int i=0; i<len-1; i++)
	{
		int move = (bitseq >> (len-i-1)) & 0x1L;
		if (0 == move)
		{
			front <<= 1;
		}
		else
		{
			front = find_leader(front);
			if (front < 0) return front;  // overflow
		}
		// printf("bit %d move=%d front=%ld\n", i, move, front);
	}
	// printf("------\n");
	return front;
}

// Do the minkowski run-length encoding trick
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
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s <max-order> <dyad-base>\n", argv[0]);
		exit(1);
	}
	int maxord = atoi(argv[1]);
	int dyad = atoi(argv[2]);

	long nmax = 1UL << maxord;
	malloc_gold(nmax);

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
		long idx = bitseq_to_idx(bits, dlen);
		if (idx < 0) { overflo++; continue; }
		if (nmax <= idx) { overflo++; continue; }

		double gold = find_gold(idx);

		double cf = bitseq_to_cf(bits, dlen);
		printf("%d	%d	%d	%ld	%g	%g	%g\n", i, bits, 1<<dlen, idx, x, gold, cf);
	}

	printf("#\n# Num overflows = %d\n#\n", overflo);
}

/* --------------------------- END OF LIFE ------------------------- */
