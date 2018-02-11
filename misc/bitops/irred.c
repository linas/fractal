/*
 * irred.c
 *
 * Find and verify irreducible golden polynomials.
 * February 2018
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>


/* Return polynomial for bitstr */
double beta(int n, double x)
{
	double acc = 0.0;
	double xn = 1.0;
	int bitstr = 2*n+1;
	while (bitstr)
	{
		if (bitstr%2 == 1) acc += xn;
		xn *= x;
		bitstr >>= 1;
	}
	return xn - acc;
}

/* Return length of bitstr, length in bits */
int len(int n)
{
	int len=0;
	int bitstr = 2*n+1;
	while (bitstr) { len++; bitstr >>= 1; }
	return len;
}

int main(int argc, char* argv[])
{
	for (int n=0; n<20; n ++)
	{
		printf("%d	%d	%g\n", n, len(n), beta(n, 1));
	}
}
