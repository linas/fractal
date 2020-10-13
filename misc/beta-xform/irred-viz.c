/*
 * irred-viz.c
 *
 * Find integer sequence for the golden polynomials.
 * ... and then 2D plot
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

void print_bitstr(int len, double gold)
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

/** Helper array, needed for finding gold midpoints */
double* zero = NULL;
void setup_gold(int nmax)
{
	zero = (double*) malloc(nmax*sizeof(double));
}

/**
 * Find the single, unique real zero of the n'th golden polynomial.
 * If the value of `n` does not really correspond to a golden
 * polynomial, return zero.
 */
double find_gold(int n)
{
	double gold = find_zero(n, 1.0, 2.0);
	zero[n] = gold;

#if 0
	bool aok = true;
	for (int j=0; j< n; j++)
	{
		double z = beta(n, zero[j]);
		if (fabs(z) < 4.0e-14) { aok = false; break; }
	}
#endif

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

	if (ork) return gold;
	return 0.0;
}

int main(int argc, char* argv[])
{
	int nmax = (1<<12) + 1;

	setup_gold(nmax);

	int cnt = 0;
	int allcnt = 0;
	int nsum = 0;
	int plen = 0;

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
			cnt = 0;
		}
		if (0.5 < gold)
		{
			cnt ++;
			allcnt ++;
			nsum += n;

			// always the case that tord/2 <= n < tord
			double frac = ((double) n) / ((double) tord);
			frac = 2.0*(frac - 0.5); // rescale to run 0 to 1.

			printf("%d	%d %d %d %g %d	%20.18g\n", allcnt, cnt, n, tord, frac, nsum, gold);
		}
	}
}
