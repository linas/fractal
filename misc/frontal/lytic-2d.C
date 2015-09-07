
#include <math.h>

#include "brat.h"
#include "lytic.h"

static double ploto(double re_q, double im_q, int itermax, double param)
{
	double x = param;

	double complex z = re_q + I * im_q;
	double complex u = z-x;
	// double complex sm = sum_extend(x, u);
	double complex sm = frac_extend(x, u);

	// return creal(sm);
	return 0.5 + 0.5 * atan2(cimag(sm), creal(sm))/M_PI;
}

DECL_MAKE_HEIGHT(ploto);
