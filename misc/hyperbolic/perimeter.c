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
		t = question_inverse (x);

		t = question_mark (t*1000000,1000000);
		t -= x;

		printf ("%d	%f	%g	%f	%f\n", i, x, t, u, q);
	}
}
