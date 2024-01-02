/*
 * recode.C
 * Recode beta orbits to 2-adic values.
 * This is the "expander" diagram, across all beta's.
 *
 * Jan 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

double recode(double beta, double x)
{
	double tn = 0.5;
	double sum = 0.0;
	for (int i=0; i<24; i++)
	{
		if (0.5 < x)
		{
			x -= 0.5;
			sum += tn;
		}
		x *= beta;
		tn *= 0.5;
	}

	return sum;
}

/*-------------------------------------------------------------------*/
/*
 */

static double coding (double x, double y, int itermax, double param)
{
	return recode(y, x);
}

DECL_MAKE_HEIGHT(coding)
