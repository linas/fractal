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

/* Reverse the bit order of i, assuming its of length n. */
int reverse(int p, int n)
{
	int r = 0;
	for (int s=0; s<n; s++)
	{
		r <<= 1;
		if (0x1 == (p & 0x1)) r |= 0x1;
		p >>= 1;
	}
	return r;
}

/*
 * argument is interpreted as p/2^n
 */
double eval (int p, int n)
{
	int m = 1<<n;
	int r = reverse (p, n);

printf("duude p=%x  r=%x\n", p,r);
	for (int i=0; i<m; i++)
	{
		
	}
	
	return 0.0;
}

int main(int argc, char * argv[])
{

	int n = 4;
	int m = 1<<n;
	for (int p=0; p<m; p++)
	{
		double x = ((double) p) / ((double) m);
		double y = eval(p, n);
		// printf("%d	%f	%f\n", p, x, y);
	} 
}
