/*
 * partition.c
 *
 * Compute the q-series-like partition function
 *
 * Linas September 2006
 */

#include <math.h>
#include <stdio.h>

double partition (double lambda)
{
	double prod = exp (lambda / (1.0-lambda));

	double x = lambda*lambda;
	double xn = 1.0;
	
	int n;
	for (n=1; n<35; n++)
	{
		xn *= x;
		double term = 1.0 + 0.5 * (exp (-xn)-1.0);
		prod *= term;

		// printf ("n=%d  term=%16.14g  prod=%16.14g\n", n, term, prod);
	}
	return prod;
}

main ()
{
	double lambda = 0.5;
	for (lambda=0.05; lambda <1.0;lambda +=0.05)
	{
		double part = partition (lambda);
		printf ("lambda=%g  prod=%16.14g\n", lambda, part);
	}
}
