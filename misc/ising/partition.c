/*
 * partition.c
 *
 * Compute the q-series-like partition function
 *
 * Linas September 2006
 */

#include <math.h>
#include <stdio.h>

main ()
{
	double lambda = 0.5;

	double prod = exp (lambda / (1.0-lambda));

	double x = lambda*lambda;
	double xn = 1.0;
	
	int n;
	for (n=1; n<40; n++)
	{
		xn *= x;
		double term = 1.0 + 0.5 * (exp (-xn)-1.0);
		prod *= term;

		printf ("n=%d  term=%g  prod=%g\n", n, term, prod);
	}
}
