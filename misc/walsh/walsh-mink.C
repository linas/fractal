/*
 * walsh-mink.C
 *
 * Dereivative of minkowski-question-mark, integrated with 
 * the Walsh functions, mapped to unit interval.
 *
 * Linas Vepstas August 2009
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

/* Reverse the bit order of i, assuming its of length n. */
int reverse(int p, int n)
{
	int r = 0;
	for (int s=0; s<n; s++)
	{
		r <<= 1;
		if (p & 0x1) r |= 0x1;
		p >>= 1;
	}
	return r;
}

/* return 1 if number of bits in p is odd. Look only at n bits. */
int odd(int p, int n)
{
	int o = 0x0;
	for (int s=0; s<n; s++)
	{
		if (p & 0x1) o = o ^ 0x1;
		p >>= 1;
	}
	return o;
}

/*
 * Given x=i/2^n, this returns the p'th walsh function 
 * at x.  i.e. returns y = w_p(x).
 *
 */
double walsh (int i, int p, int n)
{
	int r = reverse (p, n);

	int mask = r & i;
	int mo = odd (mask, n);
	double y = mo ? 1.0: -1.0;
	// printf("duude p=%x  r=%x i=%d y=%f\n", p,r, i, y);
	return y;
}

/*
 * Return integral of p'th walsh function with question mark
 *
 */
double integral (int p, int n)
{
	static ContinuedFraction f;
	f.SetEvenize();

	int m = 1<<n;
	int r = reverse (p, n);

	double acc = 0.0;

	// printf("duude p=%x  r=%x o=%d\n", p,r, o);
	double yprev = 0.0;
	for (int i=1; i<=m; i++)
	{
		int mask = r & i;
		int mo = odd (mask, n);
		
		// double x = ((double) i) / ((double) m);
		// double y = mo ? 1.0: -1.0;

		f.SetRatio(i, m);
		double y = f.ToFarey();
		
		if (mo) acc += y-yprev;
		else acc += yprev-y;

		yprev = y;
	}
	
	return acc;
}

int main(int argc, char * argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s param\n", argv[0]);
		exit(1);
	}
	int p = atoi (argv[1]);

	int n = p;
	int m = 1<<n;

	for (p=0; p<m; p++)
	{
		double x = ((double) p) / ((double) m);
		// double y = walsh(p, p, n);
		double y = integral(p, n);
		printf("%d	%f	%f\n", p, x, y);
	} 
}
