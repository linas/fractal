/*
 * Draw lines
 *
 * October 2016
 */

#include <math.h>
#include <stdio.h>
#include <time.h>

#include "brat.h"


static double ploto(double re_q, double im_q, int itermax, double param)
{
	double pixwidth = 800.0;
	double locwidth = 0.5 / pixwidth;
	double locheigh = 0.95 / pixwidth;

	double theta = M_PI * im_q;
	double x = re_q;
	double y = sin(0.5*theta);

#ifdef ALMOST
	x = pow (x, 2.5);
	double loc = param;
	if (0.9*loc < x*y and x*y < 1.1*loc) return 1.0;
	// if (loc/y - locwidth < x and x < loc/y + locwidth) return 1.0;
	// if (loc/x - locheigh < y and y < loc/x + locheigh) return 1.0;
#endif

	// Writing z = r exp(i theta) then ...
	// this curve corresponds to sqrt(r) = sine(theta/2)
	// double crv = exp(- 5.8*x);
	double crv = exp(- 5.545*x);
	if (crv - locheigh < y and y < crv + locheigh) return 1.0;

	// this curve corresponds to r/768 = sine(theta/2)
	// 768
	crv = exp(- 11.0903*(x-0.59906));
	if (crv - locheigh < y and y < crv + locheigh) return 1.0;

	// this curve corresponds to r/1200 = sine(theta/2)
	// 1200
	crv = exp(- 11.0903*(x-0.639301));
	if (crv - locheigh < y and y < crv + locheigh) return 1.0;

	// this curve corresponds to r/1500 = sine(theta/2)
	// 1500
	crv = exp(- 11.0903*(x-0.659421));
	if (crv - locheigh < y and y < crv + locheigh) return 1.0;

	// this curve corresponds to r/1948 = sine(theta/2)
	// 1948  (or 2000??)
	crv = exp(- 11.0903*(x-0.682986));
	if (crv - locheigh < y and y < crv + locheigh) return 1.0;

	crv = exp(- 11.09*(x-param));
// printf ("duuude x=%f xp=%f crv =%f\n", x, x-param, crv);
	if (crv - locheigh < y and y < crv + locheigh) return 1.0;


	return 0.0;
}


DECL_MAKE_HEIGHT(ploto);
