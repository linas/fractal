/*
 * irred-gold.c
 *
 * Given an integer index, return the corresponding "golden root" for
 * that matching beta polynomial.
 *
 * February 2018, October 2020, December 2023
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


/* Return the beta value corresponding to the n'th golden polynomial.
 * It is be constructed from the bit string of (2n+1). Construction
 * is the mid-point construction: repeated iteration of the midpoint
 * 1/2 with this beta will (re-)generate the same bitstring, until
 * returning to the midpoint. The bits are just whether the orbit went
 * left or right of midpoint. The length of the orbit will be log_2(2n+1).
 */
double gold(unsigned long n, double x)
{
	double acc = 0.0;
	double xn = 1.0;
	unsigned long bitstr = 2*n+1;
	while (bitstr)
	{
		if (bitstr%2 == 1) acc += xn;
		xn *= x;
		bitstr >>= 1;
	}
// printf("duuude n=%d x=%20.16g beta=\n", n, x, xn-acc);
	return xn - acc;
}

/* Use midpoint bisection to find the single, unique
 * positive real zero of the n'th golden polynomial.
 */
double find_zero(unsigned long n, double lo, double hi)
{
	double mid = 0.5 * (lo+hi);
	if (1.0e-15 > hi-lo) return mid;
	double fmid = gold(n, mid);
	if (0.0 < fmid) return find_zero(n, lo, mid);
	return find_zero(n, mid, hi);
}

/** Helper array, needed for finding gold midpoints. */
double* zero = NULL;
void malloc_gold(int nmax)
{
	zero = (double*) malloc((nmax+1)*sizeof(double));
	for (int i=0; i<=nmax; i++) zero[i] = -1.0;

	// Allow valid memref to zero[-1] denoting beta=1.0
	zero++;
	zero[-1] = 1.0;
	zero[0] = 2.0;
}

/** Fill up array of zero candidates. Optionally needed. */
void fill_gold(long n)
{
	// Go top down, assuming lower ranks already filled.
	for (int i=n; 0<i; i--)
	{
		if (zero[i] > -0.5) break;
		zero[i] = find_zero(i, 1.0, 2.0);
	}
}

// Return true if the polynomial root is properly bracketed for
// the index specifying that polynomial.
long zero_bracket_factor(long n, double gold)
{
	// Its valid only if it is in the middle.
#define EPS 2.0e-15
	// printf("---------\ncheck bracketing for gold=%20.16g at n=%d\n", gold, n);
	bool ork = true;
	long nhl = n;
	long nh = nhl >> 1;
	while (nh)
	{
		//printf("walk to n=%ld nhl=%ld nh=%ld znh=%g go=%g comp=%d bad=%d\n",
		//       n, nhl, nh, zero[nh], gold, 0 == nhl%2, zero[nh] <= gold);
		if (0 == nhl%2 && zero[nh] < gold+EPS) {ork = false; break; }
		nhl = nh;
		nh >>= 1;
	}
	// printf("Bracket says ork=%d\n", ork);

	if (ork) return -1;
	return nh;
}

/**
 * Return the single, unique real zero of the n'th bitstring polynomial.
 * Recall most values of `n` do NOT correspond to a golden polynomial.
 */
double find_poly_zero(long n)
{
	fill_gold(n);
	return zero[n];
}

/**
 * Return the single, unique real zero of the n'th golden polynomial.
 * If the value of `n` does not correspond to a golden polynomial,
 * return zero.
 */
double find_gold(long n)
{
	if (-1 == n) return 1.0;
	double gold = find_poly_zero(n);
	if (-1 == zero_bracket_factor(n, gold)) return gold;
	return 0.0;
}

// =================================================================

// Perform standard iteration
double tee(double beta, double x)
{
	if (x<=0.5) return beta*x;
	return beta*(x-0.5);
}

// print the bit-string that defines the orbit
void print_orbit(int len, double gold)
{
	double mid = 0.5*gold;
	while (1 < len)
	{
		if (mid < 0.5) printf(" 0");
		else printf (" 1");

		if (0.5 <= mid) mid -= 0.5;
		mid *= gold;
		len --;
	}
	printf("\n");
}

// =================================================================

/** Return order of beta polynomial for index n.
 * This is just the length of binary bitstring for 2n+1.
 */
int order(unsigned long n)
{
	int len=0;
	unsigned long bitstr = 2*n+1;
	while (bitstr) { len++; bitstr >>= 1; }
	return len;
}

/**
 * Print index as golden bitsequence.
 * Actually print order, then index, then bitseq
 */
void prt_bitstr(unsigned long n, const char* pfx, const char* sfx)
{
	int ord = order(n);
	printf("%s|%d| %ld={", pfx, ord, n);
	unsigned long bitstr = 2*n+1;
	for (int i=0; i<ord; i++)
	{
		printf("%ld", 0x1 & bitstr>>(ord-i-1));
	}
	printf("}%s", sfx);
}

// =================================================================
