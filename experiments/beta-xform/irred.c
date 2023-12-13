/*
 * irred.c
 *
 * Find and verify irreducible golden polynomials.
 * February 2018
 *
 * OBSOLETE! The code below has been copied to irred-fraction.c which
 * cleans it up, adds many additional goodies and tools and functions.
 * Don't work with this code.
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
	for (int k=0; k< nmax; k++) zero[k] = -1.0;
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

/**
 * Find the zeros, then arrange them in sequential order.
 */
int bubble_sort(int nmax)
{
	// Find all of them.
	bool bad[nmax];
	for (int n=0; n<nmax; n ++)
	{
		double gold = find_gold(n);
		bad[n] = false;
		if (gold < 0.5) bad[n] = true;
	}

	// Discard the bad ones.
	int cnt = 0;

	// XXX Also discard the early ones.
	// for (int i=0; i< nmax; i++)
	for (int i=nmax/2; i< nmax; i++)
	{
		zero[cnt] = zero[i];
		if (!bad[i]) cnt++;
	}

	// Bubble sort them.
	for (int i=0; i< cnt; i++)
	{
		for (int j=i+1; j<cnt; j++)
		{
			if (zero[j] < zero[i])
			{
				double tmp = zero[j];
				zero[j] = zero[i];
				zero[i] = tmp;
			}
		}
	}

	return cnt;
}

#ifdef BIN_COUNT
int main(int argc, char* argv[])
{
	// int nmax = (1<<28) + 3;
	int nmax = (1<<18) + 3;
	// int nmax = (1<<20) + 3;

	setup_gold(nmax);

#define NPTS 1303
	int npts = NPTS;
	double bincnt[npts+1];
	for (int i=0; i<npts; i++)
	{
		bincnt[i] = 0.0;
	}

	int cnt = 0;
	for (int n=0; n<nmax; n ++)
	{
		double gold = find_gold(n);

		// Bin-count.
		if (gold < 0.5) continue;
		cnt ++;
		int nbin = (gold - 1.0) * npts;
		bincnt[nbin] += 1.0;
	}

	double norm = ((double) npts) / ((double) cnt);
	for (int i=0; i<npts; i++)
	{
		double x = ((double) i + 0.5) / ((double) npts);
		x += 1.0;
		double y = norm * bincnt[i];
		printf("%d	%g	%g\n", i, x, y);
	}
}

#endif // BIN_COUNT

int main(int argc, char* argv[])
{
	// int nmax = (1<<28) + 3;
	// int nmax = (1<<18) + 3;
	// int nmax = (1<<20) + 3;
	// int nmax = (1<<8) + 3;
	// int nmax = (1<<12) + 1;
	int nmax = (1<<20) + 1;

	setup_gold(nmax);

// find_gold(16);
// printf("---------\ngold=%20.16g\n", zero[16]);
// exit(0);

#define PRINT_STUFF
#ifdef PRINT_STUFF
	int cnt = 0;
	int lyn = 1;
	int plen = 0;
	for (int n=0; n<nmax; n ++)
	{
		double gold = find_gold(n);

		// printf("---------\ngold=%g\n", gold);
		if (plen != len(n))
		{
			double rat = ((double) cnt) / ((double) lyn);
			printf("# total for len=%d is %d	ratio=%g\n", plen, cnt, rat);
			plen = len(n);
			lyn = cnt;
			cnt = 0;
		}
		if (0.5 < gold) { cnt++; }

		if (n < 128)
		// if (n < 540)
		// if (0)
		{
#if 0
			printf("%d l=%d %d %20.18g\n",
				n, len(n), cnt, gold);
#endif
			if (0.5 < gold)
			{
				printf("%d	%d	%20.18g #", n, len(n), gold);
				print_bitstr(len(n), gold);
			}
		}
	}
#endif


#ifdef SORTED
	printf("#\n# Columns: idx, gold, ratio, ratio-bump\n#\n");
	int cnt = bubble_sort(nmax);
	for (int n=0; n<cnt-1; n ++)
	{
		double gold = zero[n];
		double ratio = zero[n+1] / gold;
		double bump = 1.0 / (ratio - 1.0);
		printf("%d	%g %g	%g\n", n, gold, ratio, bump);
	}
#endif
}
