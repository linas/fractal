
/*
 * graph of charachteristic equation for Bernoulli map
 *
 * Linas Vepstas December 2004 
 */

#include <math.h>
#include <stdio.h>


long double 
trace (long double lam)
{
	lam *= 0.5L;

	long double acc = 0.0L;
	long double lk = lam;
	long double tk = 0.5L;
	int k;
	for (k=1; k< 500; k++)
	{
		long double term = 1.0L / ((long double) k);
		term *= lk;
		term /= 1.0L - tk;
		
		acc += term;
		if (term < 1.0e-16 * acc) break;
		lk *= lam;
		tk *= 0.5L;
	}
	acc = -acc;

	return acc;
}
 

int
main ()
{
	int i;

	int imax = 523;
	
	for (i=1; i< imax; i++)
	{
		long double x = i / ((long double) imax);

		long double tr = trace (x);
		tr = expl (tr);

		printf ("%d	%Lg 	%Lg\n", i,x, tr);
	}

	return 1;
}
