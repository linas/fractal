
/*
 * Stern-Brocot tree, new algo
 * This is the same old question mark, just a different API.
 *
 * Linas Vepstas September 2015
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Computes the Stern-Brocot tree.
 * It takes as input a dyadic fraction pos / 2^level
 * assuming 0 <= pos <= 2^level.
 * It returns the corresponding p/q value in p and q.
 *
 * In practical terms, given an input x = pos / 2^level
 * this returns the inverse question-mark, given by
 * (as per usual) ?(p/q) = pos / 2^level.
 */
void stern_brocot_tree(unsigned long pos, int level,
                       unsigned long *p, unsigned long *q)
{
	if (0 == level)
	{
		if (0 == pos)
		{
			*p = 0; *q = 1; return;
		}
		if (1 == pos)
		{
			*p = 1; *q = 1; return;
		}
	}
	if (0 == pos%2)
		return stern_brocot_tree(pos>>1, level-1, p, q);

	unsigned long pl, ql, pr, qr;
	stern_brocot_tree(pos>>1, level-1, &pl, &ql);
	stern_brocot_tree((pos>>1) + 1, level-1, &pr, &qr);
	
	*p = pl + pr;
	*q = ql + qr;
}

/**
 * Same as above, but with 128 bits
 */
void stern_brocot_tree128(unsigned __int128 pos, int level,
                          unsigned __int128* p, unsigned __int128* q)
{
	if (0 == level)
	{
		if (0 == pos)
		{
			*p = 0; *q = 1; return;
		}
		if (1 == pos)
		{
			*p = 1; *q = 1; return;
		}
	}
	if (0 == pos%2)
		return stern_brocot_tree128(pos>>1, level-1, p, q);

	unsigned __int128 pl, ql, pr, qr;
	stern_brocot_tree128(pos>>1, level-1, &pl, &ql);
	stern_brocot_tree128((pos>>1) + 1, level-1, &pr, &qr);
	
	*p = pl + pr;
	*q = ql + qr;
}

int pr128(unsigned __int128 val)
{
	if (val <= ULONG_MAX)
	{
		unsigned long u64 = val;
		return printf("%lu", u64);
	}

/*       ULONG_MAX 18446744073709551615UL    */
#define P10_UINT64 10000000000000000000UL   /* 19 zeroes */
	unsigned __int128 lead = val / P10_UINT64;
	unsigned __int128 trail = val % P10_UINT64;
	unsigned long leading = lead;
	unsigned long trailing = trail;
	return printf("%lu%lu", leading, trailing);
}

#ifdef TEST
int main(int argc, char *argv[])
{

	int level = atoi(argv[1]);

	unsigned long p;
	unsigned long q;
	double norm = 1.0 / (double)(1<<level);
	for (int i=0; i<= (1<<level); i++)
	{
		stern_brocot_tree(i, level, p, q);
		double x = i * norm;
		double y = ((double) p) / (double) q;
		printf("%d	%lu	%lu	%g	%g\n", i, p, q, x, y);
	}
}
#endif
