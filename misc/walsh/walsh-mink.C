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

double question(int p, int q)
{
	static ContinuedFraction f;
	f.SetEvenize();

	f.SetRatio(p, q);
	double y = f.ToFarey();
	if (p==q) y = 1.0;

	return y;
}

double triangle (double x)
{
	x -= floor(x);
	if (x < 0.5) return x;
	return 1.0-x;
}

double takagi(int p, int q)
{
	double w = 0.6;
	double x = ((double) p) / ((double) q);
	double acc = 0.0;
	double tw = 1.0;
	double wn = 1.0;
	while(1)
	{
		acc += wn * triangle(tw*x);
		wn *= w;
		tw *= 2.0;
		if (wn < 1.0e-12) break;
	}
	return acc;
}

double line(int p, int q)
{
	double x = ((double) p) / ((double) q);
	return x;
}


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
double integral (int p, int n, double (*f)(int, int))
{
	int m = 1<<n;
	int r = reverse(p, n);
	int o = odd(p, n);

	double acc = 0.0;

	// printf("duude p=%x  r=%x o=%d\n", p,r, o);
	double yprev = f(0,m);
	for (int i=1; i<=m; i++)
	{
		int mask = r & (i-1);
		int mo = o ^ odd (mask, n);  // I think this is right ... 
		
		// double x = ((double) i) / ((double) m);
		// double w = mo ? 1.0: -1.0;

		double y = f(i, m);

		if (mo) acc += y-yprev;
		else acc += yprev-y;

		// printf("%f	%f	%f %f\n", x,w,y,acc);
		yprev = y;
	}
	
	return -acc;
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

	// integral (3, 4);

#define UNIT_INTERVAL 1
#if UNIT_INTERVAL
	// double (*f)(int, int) = question;
	// double (*f)(int, int) = line;
	double (*f)(int, int) = takagi;

	double acc = 0.0;
	for (p=1; p<m; p++)
	{
		double x = ((double) p) / ((double) m);
		// double y = walsh(p, p, n);
		// r is the r'th Walsh function.
		int r = reverse (p, n);
		double y = integral(r, n, f);
		acc += y;
		double fx = f(p,m);
		printf("%d	%f	%f	%f	%f\n", p, x, y, acc, fx);
	} 
#endif

// #define SERIES 1
#ifdef SERIES
	double acc = 0.0;
	for (p=1; p<m; p++)
	{
		double x = p;
		double y = integral(p, n);
		acc += y;
		printf("%d	%f	%f	%f\n", p, x, y, acc);
	} 
#endif
}
