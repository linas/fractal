/* 
 * perimeter.c
 *
 * Map tree to perimeter of the circle.  Used to sanity check
 * various results.
 *
 * Linas Vepstas May 2007
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "flt.h"
#include "question.h"

/**
 * dyadic_to_stern_brocot -- convert real number to stern-brocot
 *
 * This function maps the unit interval, expressed as a dyadic
 * tree, to the full stern brocot tree.  The input is 0<x<1.
 * Perform binary expansion of x. Use L,R for binary digits 0,1. 
 * Concatenate the L and R to form a matrix. Padd on the right with
 * zeros. The resulting fractinal linear transform will map the 
 * entire complex plane to b/d. This function returns b/d.
 *
 * Basically, this function maps the dyadic tree to the 
 * Stern-Brocot tree.
 */
double dyadic_to_stern_brocot (double x)
{
	if (0.0>x) x -= (int) x - 1;
	if (1.0<x) x -= (int) x;

	mobius_t ell = mobius_set (1,0,1,1);
	mobius_t are = mobius_set (1,1,0,1);

	mobius_t acc = mobius_set (1,0,0,1);
	int i;
	for (i=0; i<45; i++)
	{
		if (0.5 <= x)
		{
			acc = mobius_mul (acc, are);
			x -= 0.5;
		}
		else
		{
			acc = mobius_mul (acc, ell);
		}
		x *= 2.0;
	}
	return acc.b.re / acc.d.re;
}

/**
 * qinverse - return the inverse of the question mark function
 *
 * This impments a rapid algorithm to compute the inverse of 
 * the question mark function.
 */
double qinverse (double x)
{
	return dyadic_to_stern_brocot (0.5*x);
}

main()
{
	int i;

	int deno=32;
	for (i=0; i<deno;i++)
	{
		double x = ((double) i)/((double) deno);

		double t = dyadic_to_stern_brocot(x);

		double u = atan2 (-2.0*t, t*t-1.0);
		// double u = atan2 (-4.0*t, t*t-4.0);
		// double u = atan2 (-t, t*t-0.25);

		// if (0.0> u) u += 2.0*M_PI;
		// u /= 2.0*M_PI;
		u /= M_PI;
		u += 1.0;
		
		// double q = question_mark (i,deno);
		// double q = question_inverse (2.0*x);
		// q = 0.5 + 0.5*q;
		// x = 0.5 + 0.5*x;

#if 0
		double h = t;
		if (h > 1.0) h = 1.0/h;
		h = question_mark (h*1000000,1000000);
		h *= 0.5; 
		if (t>1.0)
		{
			h = 1.0-h;
		}
#endif
	
		double q=0;
#if 0
		if (x<0.5)
		{
			q = question_inverse (2.0*x);
		}
		else
		{
			q = 1.0/question_inverse (2.0*(1.0-x));
		}
		q = question_inverse (x);
#endif
		t = qinverse (x);

		t = question_mark (t*1000000,1000000);
		t -= x;

		printf ("%d	%f	%g	%f	%f\n", i, x, t, u, q);
	}
}
