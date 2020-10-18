/*
 * irred-fraction.c
 *
 * Find integer sequence for the golden polynomials.
 * Relate it to the continued-fraction representation.
 * See irred.c for additional utilities.
 *
 * February 2018, October 2020
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


/* Return the n'th golden polynomial. It can be constructed from
 * the bit string of (2n+1).
 */
double beta(unsigned long n, double x)
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
	double fmid = beta(n, mid);
	if (0.0 < fmid) return find_zero(n, lo, mid);
	return find_zero(n, mid, hi);
}

/* Return length of bitstr, length in bits */
int len(unsigned long n)
{
	int len=0;
	unsigned long bitstr = 2*n+1;
	while (bitstr) { len++; bitstr >>= 1; }
	return len;
}

/** Helper array, needed for finding gold midpoints */
double* zero = NULL;
void setup_gold(int nmax)
{
	zero = (double*) malloc(nmax*sizeof(double));
	for (int i=0; i< nmax; i++) zero[i] = -1.0;
}

// Return true if this is a valid polynomial
bool zero_is_bracketed(int n, double gold)
{
	// Its valid only if it is in the middle.
#define EPS 2.0e-15
	// printf("---------\ngold=%g\n", gold);
	bool ork = true;
	int nhl = n;
	int nh = nhl >> 1;
	while (nh)
	{
		// printf("duuude n=%d nhl=%d nh=%d znh=%g go=%g comp=%d bad=%d\n",
		// n, nhl, nh, zero[nh], gold, 0 == nhl%2, zero[nh] <= gold);
		if (0 == nhl%2 && zero[nh] < gold+EPS) {ork = false; break; }
		nhl = nh;
		nh >>= 1;
	}

	return ork;
}

/**
 * Find the single, unique real zero of the n'th golden polynomial.
 * If the value of `n` does not really correspond to a golden
 * polynomial, return zero.
 */
double find_gold(int n)
{
	if (-1 == n) return 2.0;

	double gold = zero[n];
	if (gold < -0.5)
	{
		gold = find_zero(n, 1.0, 2.0);
		zero[n] = gold;
	}

	if (zero_is_bracketed(n, gold)) return gold;
	return 0.0;
}

/*
 * Given the continue-fraction representation, return the
 * corresponding sequence number.
 */
long sequence_from_cf(int cfrac[], int len)
{
	if (1 == len)
	{
		if (-1 == cfrac[0]) return -1;
		return 1 << cfrac[0];
	}

	if (3 < len) return -1; // unknown.

	long leader = sequence_from_cf(cfrac, len-1);

	long follower = 2*leader + 1;
	follower *= 1 << (cfrac[0] + cfrac[len-1]);

	return follower;
}

void print_seq(int cfrac[], int len, char* head, char* tail)
{
	printf("%s [", head);
	for (int i=0; i<len; i++) printf(" %d", cfrac[i]);
	printf("]%s", tail);
}

#define SZ 10
/*
 * Iterate on the continued fraction.
 * i.e. generate sequences
 * maxdepth == number of doubling steps
 * maxlength == max length of fraction.
 * maxn == cutoff fir highest known n
 */
void iterate_cf(int cfrac[], int len, int maxdepth, int maxlength, long maxn)
{
	long seq = sequence_from_cf(cfrac, len);
	if (seq >= maxn) return;
	double gold = find_gold(seq);
	printf("seq = %ld gold=%g ", seq, gold);
	print_seq(cfrac, len, "", "");

	// Validate bracketing.
	if (1 < len)
	{
		print_seq(cfrac, len-1, " left=", "");

		long left = sequence_from_cf(cfrac, len-1);

		// Right bracket is tricky.
		int rfrac[SZ];
		for (int i=0; i<len; i++) rfrac[i] = cfrac[i];
		rfrac[len-2]--;
		print_seq(rfrac, len-1, " right=", "");

		long right = sequence_from_cf(rfrac, len-1);

		if (left <= right)
			printf("FAIL bracket order! left=%ld right=%ld\n", left, right);

		if (seq <= left)
			printf("FAIL left bracket! left=%ld seq=%ld\n", left, seq);

		if (seq <= right)
			printf("FAIL right bracket! right=%ld seq=%ld\n", right, seq);

		double lg = find_gold(left);
		double rg = find_gold(right);

		if (gold < lg) printf("fail left gold! left=%g gold=%g\n", lg, gold);
		if (rg < gold) printf("fail right gold! gold=%g right=%g\n", gold, rg);
	}

	printf("\n");

	// Iterate length-wise first
	if (len < maxlength)
	{
		int bfrac[SZ];
		for (int i=0; i<len; i++) bfrac[i] = cfrac[i];
		bfrac[len] = 0;
		iterate_cf(bfrac, len+1, maxdepth, maxlength, maxn);
	}

	// Iterate depthwise second.
	if (cfrac[len-1] < maxdepth)
	{
		int bfrac[SZ];
		for (int i=0; i<len; i++) bfrac[i] = cfrac[i];
		bfrac[len-1] ++;
		iterate_cf(bfrac, len, maxdepth, maxlength, maxn);
	}
}

int main(int argc, char* argv[])
{
	int nmax = (1<<10) + 1;

	setup_gold(nmax);

	for (int n=0; n<nmax; n ++)
		find_gold(n);

	int cfrac[SZ];
	cfrac[0] = 0;

	iterate_cf(cfrac, 1, 2, 2, nmax);
}
