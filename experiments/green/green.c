
/*
 * green.c 
 *
 * Green functions for the Bernoulli process
 *
 * Linas Vepstas December 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "binomial.h"
#include "bernoulli.h"
#include "gaussian.h"

double do_sum (Gaussian *g, double x)
{
	int i;

	double acc = 0.0;
	double tn = 2.0;
	for (i=1; i<16; i++)
	{
		double right = gaussian_eval (g, i-1, 1.0);
		right -= gaussian_eval (g, i-1, 0.0);
		right /= factorial (i);

		double left = bernoulli_poly (i, x);
		acc += right * left * tn;
		tn *= 2.0;
	}

	return acc;
}

main (int argc, char * argv[])
{
	double x;

	double sigma = 0.15;
	double mu = 0.4;

	sigma = atof (argv[1]);

	Gaussian *g = gaussian_new (mu, sigma);

	for (x=0.0; x<=1.0; x+= 0.01)
	{
		double y = do_sum (g, x);
		printf ("%g	%g\n", x, y);
	}
}

