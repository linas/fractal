/*
 * selfie-util.c
 * Misc utilities for selfie.c
 *
 * December 2023
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

/* ================================================================= */

/* Return greatest common divisor of p and q */
unsigned long gcd(unsigned long p, unsigned long q)
{
	/* Euclid's algorithm for obtaining the gcd */
	while (0 != q)
	{
		unsigned long t = p % q;
		p = q;
		q = t;
	}
	return p;
}

/* ================================================================= */
// Convert dyadic fraction to canonical tree node number.
// The canonical tree numbering is 1 at the root, 2,3 in first row,
// 4,5,6,7 in next row, etc.
// 
// A dyadic fraction is in the form p/q with p=2m+1 and q=2^n
// The conversion is then m + 2^(n-1) = (p+q-1)/2
unsigned long dyafrac_to_tree_idx(unsigned long p, unsigned long q)
{
	unsigned long gcf = gcd(p, q);
	p /= gcf;
	q /= gcf;
	return (p + q) >> 1;
}

void tree_idx_to_dyafrac(unsigned long idx, unsigned long* p, unsigned long* q)
{
	int len = bitlen(idx);
	unsigned long deno = 1UL << len;
	unsigned long em = idx - (deno>>1);
	unsigned numo = 2*em + 1UL;
	*p = numo;
	*q = deno;
}

/* Return canonical dyadic associated with canonical integer
 * in the canonical tree numbering.The map is
 *   1 |--> 1/2
 *   2 |--> 1/4
 *   3 |--> 3/4
 *   4 |--> 1/8
 *   5 |--> 3/8
 * etc.
 */
double canonical_dyadic(unsigned long idx)
{
	unsigned long p=0, q=0;
	tree_idx_to_dyafrac(idx, &p, &q);
	return ((double) p) / ((double) q);

#ifdef ALT_IMPLEMENTATION
	int len = bitlen(n);
	unsigned long base = 1UL<<len;
	unsigned long nb = n - base/2;
	unsigned long nd = 2*nb+1;
	double dya = ((double) nd) / ((double) base);
	return dya;
#endif
}

/* --------------------------- END OF LIFE ------------------------- */
