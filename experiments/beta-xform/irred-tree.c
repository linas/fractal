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

void print_bitseq(long bitseq, int len, char* pre, char* suf)
{
	printf("%s", pre);
	for (int i=0; i<len; i++)
	{
		int bit = (bitseq >> (len-i-1)) & 0x1L;
		printf("%d", bit);
	}
	printf (" /%d", len);
	printf("%s", suf);
}

// Get the index is the left side of the bracket
long get_left_idx(long idx)
{
	while (0 == idx%2)
		idx >>= 1;

	// And once more
	idx = (idx - 1L) / 2;
	return idx;
}

// Given an index, return the matching tree-walk to get to it.
// Inverse of bitseq_to_idx
long idx_to_bitseq(long idx, int* leng)
{
	*leng = 0;
	int len = 0;
	long bitseq = 0;
	// printf("enter idx=%ld\n", idx);
	while (0 < idx)
	{
		bitseq |= 1UL << len;
		len++;
		while (0 == idx%2 && is_valid_index(idx/2UL))
		{
			idx >>= 1;
			len++;
		}

		long left = get_left_idx(idx);
		// printf("%ld |=> leader idx=%ld len=%d\n", left, idx, len);
		idx = left;
	}
	*leng = len;
	return bitseq;
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
#ifdef SPOT_CHECK
	int idx = atoi(argv[1]);
	malloc_gold(idx+1);

	int len = 0;
	long bits = idx_to_bitseq(idx, &len);
	print_bitseq(bits, len, " # bits=(", ")\n");

	long ridx = bitseq_to_idx(bits, len);
	printf("reconstruct %ld\n", ridx);
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
		if (false == is_valid_index(idx)) continue;

		int len = 0;
		long bits = idx_to_bitseq(idx, &len);
		// print_bitseq(bits, len, " # bits=(", ")\n");

		long ridx = bitseq_to_idx(bits, len);
		// printf("reconstruct %ld\n", ridx);

		if (ridx != idx)
		{
			printf("Error: bad indexing at idx=%ld rec=%ld", idx, ridx);
			print_bitseq(bits, len, " # bits=(", ")\n");
		}
	}
#endif

#define PRINT_DYADIC_TO_BETA_MAP
#ifdef PRINT_DYADIC_TO_BETA_MAP
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

		// And now, a shifted version
		int shbit = bits;
		long loidx = bitseq_to_idx(shbit, dlen+1);
		if (loidx < 0) { overflo++; continue; }
		if (nmax <= loidx) { overflo++; continue; }
		double logold = find_gold(loidx);

		shbit = bits | 1<<dlen;
		long hiidx = bitseq_to_idx(shbit, dlen+1);
		double higold = find_gold(hiidx);

		double cf = bitseq_to_cf(bits, dlen);
		printf("%d	%d	%d	%ld	%g	%g	%g", i, bits, 1<<dlen, idx, x, gold, cf);
		printf("	%g", logold);
		printf("	%g", higold);

		print_bitseq(bits, dlen, " # bits=(", ")\n");
	}

	printf("#\n# Num overflows = %d\n#\n", overflo);
#endif
}

/* --------------------------- END OF LIFE ------------------------- */
