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

double delta(double x)
{
	if (0.0>x) x = 0.0;
	if (0.25>x) return -1.0/question_inverse (4.0*x);
	if (0.5>x) return -question_inverse (2.0-4.0*x);
	if (0.75>x) return question_inverse (4.0*x-2.0);
	return 1.0/question_inverse (4.0-4.0*x);
}

main()
{
	int i;

	int deno=443;
	for (i=0; i<deno;i++)
	{
		double x = ((double) i)/((double) deno);

		// double t = dyadic_to_stern_brocot(x);
		double t = delta(x);

		double theta = atan2 (-2.0*t, t*t-1.0);
		double phi = atan2 (-4.0*t, t*t-4.0);
		// double u = atan2 (-t, t*t-0.25);

		if (0.0> theta) theta += 2.0*M_PI;
		theta /= 2.0*M_PI;
		
		if (0.0> phi) phi += 2.0*M_PI;
		phi /= 2.0*M_PI;
		
		// double q = question_mark (i,deno);
		double q = question_inverse (x);

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
	
		printf ("%d	%f	%f	%f	%f	%f\n", i, x, t, theta, phi, q);
	}
}
