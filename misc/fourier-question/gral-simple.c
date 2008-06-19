/*
 * Compute matrix elements of fourier of minkowski question mark 
 * by means of brute-force integration. This seems to be a viable
 * quick-n-dirty way of getting these.
 *
 * Linas June 2008
 */

#define M_PIl    3.1415926535897932384626433832795029L  /* pi */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void make_elt(int m, int n, long double *pre, long double *pim)
{
	long double start, end, th, off, g;
	long double re = 0.0L;
	long double im = 0.0L;
	int i, j;


	*pre = re;
	*pim = im;
}

void fill_matrix(int sz)
{
	int i, j;

	for (i=-sz; i<=sz; i++)
	{
		for (j=-sz; j<=sz; j++)
		{
			long double re = 0.0L;
			long double im = 0.0L;

			make_elt(i,j, &re, &im);
			printf("W[%d, %d] = %Lg +i %Lg\n", i, j, re, im);
		}
		printf ("\n");
	}
}

main()
{
	fill_matrix(3);
}
