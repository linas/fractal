
/*
 * galois.c
 *
 * solution to galois-theory inspired schroedinger equation
 * January 2006 -- Linas Vepstas
 */

#include <stdlib.h>

/* Return a_n given higer an's and energy */
double 
zero_ansatz (int n, double a8, double a6, double a4, double a2, double e)
{
	double a = a8 *(n*n+15*n+56);
	a += a6*(e -8*n*n -8*13*n -16*21-14-0.75);
	a += a4*(16*n*n + 16*11*n +16*28 + 5 -16*e);
	a += a2*(16*e - 36);

	a /= 16.0;
	return a;
}

main(int argc, char* argv[])
{
	int n;
#define NMAX 1000
	double ar[NMAX];

	ar[NMAX-2] = 1.0;
	ar[NMAX-4] = 0.5;
	ar[NMAX-6] = 3.0/16.0;
	ar[NMAX-8] = 1.0/16.0;
	for (n=NMAX-10; n>=0; n-=2)
	{
		ar[n] = zero_ansatz (n, ar[n+8], ar[n+6], ar[n+4], ar[n+2], e);
		printf ("%d	%g\n", n, ar[n]);
	}
}

