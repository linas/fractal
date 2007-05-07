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
 * unit_interval_to_mobius -- convertt real number to mobius xform
 *
 * Innput is real number x, assume 0<x<1.  Perform binary 
 * expansion of x. Use L,R for binary digits 0,1.  Return 
 * the result.
 */
mobius_t unit_interval_to_mobius (double x)
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

	/* assume terminated by infinite string of zeros */
	mobius_t lim = mobius_set_d (1,0,1.0e200,1);
	acc = mobius_mul (acc, lim);

	return acc;
}

double mobius_to_limit (mobius_t m)
{
	return m.b.re/m.d.re;
}

main()
{
	int i;

	int deno=32;
	for (i=0; i<deno;i++)
	{
		double x = ((double) i)/((double) deno);

		mobius_t m = unit_interval_to_mobius (x);
		double t = mobius_to_limit(m);

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

		double h = t;
		if (h > 1.0) h = 1.0/h;
		h = question_mark (h*1000000,1000000);
		h *= 0.5; 
		if (t>1.0)
		{
			h = 1.0-h;
		}
	
		double q;
		if (x<0.5)
		{
			q = question_inverse (2.0*x);
		}
		else
		{
			q = 1.0/question_inverse (2.0*(1.0-x));
		}

		//printf ("%d	%f	%f	%f	%f\n", i, x, t, u, q);
	}
}
