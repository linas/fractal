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
	double locheigh = 0.75 / pixwidth;

	double theta = M_PI * im_q;
	double x = re_q - 1.0;
	double y = sin(0.5*theta);

	double loc = 0.1;
	if (0.9*loc < x*y and x*y < 1.1*loc) return 1.0;

	if (loc/y - locwidth < x and x < loc/y + locwidth) return 1.0;
	if (loc/x - locheigh < y and y < loc/x + locheigh) return 1.0;

	return 0.0;
}


DECL_MAKE_HEIGHT(ploto);
