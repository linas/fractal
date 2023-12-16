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
	for (int i=1; i< len; i++)
	{
		int move = (bitseq >> (len-i-1)) * 0x1L;
		if (move)
		{
			front <<= 1;
		}
		else
		{
			front = find_leader(front);
			if (front < 0) return front;  // overflow
		}
	}
	return front;
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
		double x = ((double) i - 0.5) / (double) maxdy;

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
		printf("%d	%d	%d	%ld	%g	%g\n", i, bits, dlen, idx, x, gold);
	}

	printf("#\n# Num overflows = %d\n#\n", overflo);
}

/* --------------------------- END OF LIFE ------------------------- */
