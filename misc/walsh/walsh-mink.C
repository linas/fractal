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
	int m = 1<<n;
	int r = reverse (p, n);

	int mask = r & i;
	int mo = odd (mask, n);
	double y = mo ? 1.0: -1.0;
	// printf("duude p=%x  r=%x i=%d y=%f\n", p,r, i, y);
	return y;
}

/*
 * argument is interpreted as p/2^n
 */
double eval (int p, int n)
{
	int m = 1<<n;
	int r = reverse (p, n);
	int o = odd(p, n);

	// printf("duude p=%x  r=%x o=%d\n", p,r, o);
	for (int i=0; i<m; i++)
	{
		int mask = r & i;
		int mo = odd (mask, n);
		
		double x = ((double) i) / ((double) m);
		double y = mo ? 1.0: -1.0;

		printf("%d	%f	%f\n", i, x, y);
	}
	
	return 0.0;
}

int main(int argc, char * argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s param\n", argv[0]);
		exit(1);
	}
	int p = atoi (argv[1]);

	int n = 4;
	int m = 1<<n;

	eval(p, n);
#if 0
	for (int p=0; p<m; p++)
	{
		double x = ((double) p) / ((double) m);
		double y = eval(p, n);
		// printf("%d	%f	%f\n", p, x, y);
	} 
#endif
}
