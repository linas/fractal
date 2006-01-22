
/*
 * galois.c
 *
 * solution to galois-theory inspired schroedinger equation
 * January 2006 -- Linas Vepstas
 */

#include <stdio.h>
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

	if (argc != 2)
	{
		fprintf (stderr, "Usage: %s <energy>\n", argv[0]);
		exit (1); 
	}
	double e = atof (argv[1]);

	int nlast = 100;

	ar[nlast] = 0.0;
	ar[nlast-2] = 0.0;
	ar[nlast-4] = 0.0;
	ar[nlast-6] = 1.0;
	for (n=nlast-8; n>=0; n-=2)
	{
		ar[n] = zero_ansatz (n, ar[n+8], ar[n+6], ar[n+4], ar[n+2], e);
	}

	int nstart = 0;
	if (nlast%2)
	{
		nstart = 1;
	}

	/* normalize */
	double off = 1.0/ar[nstart];
	for (n=nstart; n<nlast; n+=2)
	{
		ar[n] *= off;
		// printf ("%d	%g\n", n, ar[n]);
	}

	/* wave function */
	double x;

	for (x=-10.0; x<10.0; x +=0.4)
	{
		double xn = 1.0;
		if (nlast%2)
		{
			xn = x;
		}
		double psi = 0.0;
		for (n=nstart; n<nlast; n+=2)
		{
			psi += ar[n] * xn;
			xn *= x*x;
		}
		psi *= x*x-0.25;
		printf ("%g	%g\n", x, psi);
	}
}

