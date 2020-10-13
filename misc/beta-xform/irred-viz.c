/*
 * irred-viz.c
 *
 * Find integer sequence for the golden polynomials.
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
	double gold = zero[n];
	if (gold < -0.5)
	{
		gold = find_zero(n, 1.0, 2.0);
		zero[n] = gold;
	}

	if (zero_is_bracketed(n, gold)) return gold;
	return 0.0;
}

int main(int argc, char* argv[])
{
	int nmax = (1<<21) + 1;

	setup_gold(nmax);

	int cnt = 0;
	int allcnt = 0;
	int nsum = 0;
	int plen = 0;
	int lyn = 1;

	int tord = 1<<plen;

	for (int n=0; n<nmax; n ++)
	{
		double gold = find_gold(n);

		// printf("---------\ngold=%g\n", gold);
		if (plen != len(n))
		{
			printf("# total for len=%d is %d\n", plen, cnt);
			plen = len(n);
			tord = 1<<plen;
			tord /= 2;
			lyn = cnt;
			cnt = 0;

			// advance and count lyndon number
			lyn = 0;
			for (int k=n+1; k<nmax; k++)
			{
				double g = find_gold(k);
				if (0.5 < g) lyn++;
				if (len(k) != plen) break;
			}
			printf("# next lyndon is %d\n", lyn);
		}
		if (0.5 < gold)
		{
			cnt ++;
			allcnt ++;
			nsum += n;

			// always the case that tord/2 <= n < tord
			double frac = ((double) n) / ((double) tord);
			frac = 2.0*(frac - 0.5); // rescale to run 0 to 1.

			// OK, this should work.
			double seqfrac = ((double) cnt) / ((double) lyn);

			printf("%d	%d %d %d %g %g	%d	%20.18g\n",
				allcnt, cnt, n, tord, seqfrac, frac, nsum, gold);
		}
	}
}
