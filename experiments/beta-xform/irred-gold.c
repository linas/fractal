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
double golden_poly(unsigned long n, double x)
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
	double fmid = golden_poly(n, mid);
	if (0.0 < fmid) return find_zero(n, lo, mid);
	return find_zero(n, mid, hi);
}

/** Helper array, needed for finding gold midpoints. */
static double* zero = NULL;
static int* stopper = NULL;
static long maxidx = -2;
void malloc_gold(long nmax)
{
	if (-2 != maxidx)
	{
		printf("Errror: gold array already malloced\n");
		return;
	}
	maxidx = nmax;
	zero = (double*) malloc((nmax+1)*sizeof(double));
	for (int i=0; i<=nmax; i++) zero[i] = -1.0;
	stopper = (int*) malloc((nmax+1)*sizeof(int));
	for (int i=0; i<=nmax; i++) stopper[i] = 0;

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
long theta_factor(long n, double gold)
{
	// Its valid only if it is in the middle.
#define EPS 2.0e-15

// #define DBZ(X) printf X
#define DBZ(X)
	DBZ(("---------\ncheck bracketing for gold=%20.16g at n=%ld\n", gold, n));
	bool ork = true;
	long nhl = n;
	long nh = nhl >> 1;
	while (nh)
	{
		DBZ(("walk to n=%ld nhl=%ld nh=%ld znh=%g go=%g comp=%d bad=%d\n", \
		     n, nhl, nh, zero[nh], gold, 0 == nhl%2, zero[nh] <= gold));

		// if (0 == nhl%2 && zero[nh] < gold+EPS) stopper[nh] ++;
		if (0 == nhl%2 && zero[nh] < gold+EPS) {ork = false; break; }
		nhl = nh;
		nh >>= 1;
	}
	DBZ(("Bracket says ork=%d\n", ork));

	if (ork) return -1L;
	return nh;
}

bool theta(long n, double gold)
{
	return -1L == theta_factor(n, gold);
}

// Same as above; but remember result from last time, so that
// we don't compare to any non-golden roots!
// Return true if the polynomial root is properly bracketed for
// the index specifying that polynomial.
long stopit(long n, double gold)
{
	// Its valid only if it is in the middle.
#define EPS 2.0e-15

#define DBS(X) printf X
// #define DBS(X)
	DBS(("---------\ncheck stopid for gold=%20.16g at n=%ld\n", gold, n));
	bool ork = true;
	long nhl = n;
	long nh = nhl >> 1;
	while (nh)
	{
		DBS(("walk to n=%ld nhl=%ld nh=%ld znh=%g go=%g comp=%d bad=%d\n", \
		     n, nhl, nh, zero[nh], gold, 0 == nhl%2, zero[nh] <= gold));

		if (0 == nhl%2)
		{
			if (stopper[nh]) { ork = false; break; }
			if (zero[nh] < gold+EPS) { ork = false; break; }
		}
		nhl = nh;
		nh >>= 1;
	}
	DBS(("Bracket says ork=%d\n", ork));

	if (ork) return -1L;
printf("duude %ld stopped by %ld stopper=%d\n", n, nh, stopper[nh]);
	if (stopper[nh]) stopper[nh]++;
	stopper[n]++;
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
	// Allow silent out-of-bounds
	if (maxidx <= n) return 0.0;
	if (n < -1L) return 0.0;

	if (-1L == n) return 1.0;
	double gold = find_poly_zero(n);
	if (-1L == theta_factor(n, gold)) return gold;
	return 0.0;
}

/**
 * Return true if `n` corresponds to a golden polynomial, else false.
 */
bool is_valid_index(long n)
{
	// Allow silent out-of-bounds
	if (maxidx <= n) return false;
	if (n < -1L) return false;

	if (-1L == n) return true;
	double gold = find_poly_zero(n);
	return (-1 == theta_factor(n, gold));
}

bool is_stopper(long n)
{
	double gold = find_poly_zero(n);
	return -1L != stopit(n, gold);
}

void print_stoppers(long nmax)
{
	for (long n=0; n<nmax; n++)
		printf("%ld	%d\n", n, stopper[n]);
}

// =================================================================

// Return the n'th element of the "valid index sequence".
// This is the sequence 1,2,3,4,6,7,8,10,12,13,14,15,16,24...
// The composition is_valid_index(valid_index(n)) always returns true.
// All valid indexes are included in the list.
//
long valid_index(long n)
{
	long idx = 1;
	long cnt = 0;
	while (cnt < n)
	{
		if (is_valid_index(idx)) cnt++;
		idx++;
	}
	idx --;
	// printf("counted %ld is %ld\n", n, idx);
	return idx;
}

static long* akk = NULL;
void malloc_index_cache(long maxseq)
{
	akk = (long *) malloc(maxseq * sizeof(long));
	for (long i=0; i< maxseq; i++) akk[i] = -1;
	akk[0] = 1;
	akk[1] = 1;
}

// Caching version of above, avoids recompuation.
long valid_index_cache(long n)
{
	if (NULL == akk)
	{
		printf("Error: Failed to initialize index cache\n");
		return -2;
	}

	if (akk[n] < 0)
	{
		long idx = valid_index_cache(n-1) + 1;
		long cnt = n-1;
		while (cnt < n)
		{
			if (is_valid_index(idx)) cnt++;
			idx++;
		}
		idx --;
		akk[n] = idx;
		// printf("allowed %ld is %ld\n", n, idx);
	}

	return akk[n];
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
#define DEPS 1e-14
	double mid = 0.5*gold;
	while (1 < len)
	{
		if (mid < 0.5-DEPS) printf(" 0");
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
