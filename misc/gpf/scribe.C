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
	// double locwidth = 1.2 / pixwidth;
	double locheigh = 1.2 / pixwidth;

	// if (0.25 < im_q) locheigh = 2.0 / pixwidth;
	// if (0.5 < im_q) locheigh = 4.0 / pixwidth;
	// if (0.75 < im_q) locheigh = 8.0 / pixwidth;
	locheigh = 2.0* pow(8.0, im_q) / pixwidth;

	double theta = M_PI * im_q;
	double x = re_q;
	double y = sin(0.5*theta);

	int ipix = (int) (im_q*pixwidth/ 20.0);
	int dash = 1 - ipix%2;

#ifdef ALMOST
	x = pow (x, 2.5);
	double loc = param;
	if (0.9*loc < x*y and x*y < 1.1*loc) return 1.0;
	// if (loc/y - locwidth < x and x < loc/y + locwidth) return 1.0;
	// if (loc/x - locheigh < y and y < loc/x + locheigh) return 1.0;
#endif

	// Writing z = r exp(i theta) then ...
	// this curve corresponds to sqrt(r) = 1/sine(theta/2)
	// double crv = exp(- 5.8*x);
	double crv = exp(- 5.545*x);
	if (dash && crv - locheigh < y and y < crv + locheigh) return 1.0;

	// this curve corresponds to 768/r = sine(theta/2)
	// 768
	crv = exp(- 11.0903*(x-0.59906));
	if (dash && crv - locheigh < y and y < crv + locheigh) return 1.0;

	// this curve corresponds to 1200/r = sine(theta/2)
	// 1200
	crv = exp(- 11.0903*(x-0.639301));
	if (dash && crv - locheigh < y and y < crv + locheigh) return 1.0;

	// this curve corresponds to 1500/r = sine(theta/2)
	// 1500
	crv = exp(- 11.0903*(x-0.659421));
	if (dash && crv - locheigh < y and y < crv + locheigh) return 1.0;

	// this curve corresponds to 1948/r = sine(theta/2)
	// 1948  (or 2000??)
	crv = exp(- 11.0903*(x-0.682986));
	if (dash && crv - locheigh < y and y < crv + locheigh) return 1.0;

/*
	crv = exp(- 11.09*(x-param));
// printf ("duuude x=%f xp=%f crv =%f\n", x, x-param, crv);
	if (crv - locheigh < y and y < crv + locheigh) return 1.0;
*/


	return 0.0;
}


DECL_MAKE_HEIGHT(ploto);
