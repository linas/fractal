
/*
 * green.c 
 *
 * Green functions for the Bernoulli process
 *
 * Linas Vepstas December 2005
 */

#include <math.h>
#include <stdio.h>
#include "bernoulli.h"

main ()
{
	double x;

	for (x=0.0; x<=1.0; x+= 0.1)
	{
		double b0 = bernoulli_poly (0, x);
		double b1 = bernoulli_poly (1, x);
		double b2 = bernoulli_poly (2, x);
		double b3 = bernoulli_poly (3, x);
		double b4 = bernoulli_poly (4, x);

		printf ("%g	%g	%g	%g	%g	%g\n", x, b0, b1, b2, b3, b4);
	}
}

