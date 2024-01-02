/*
 * recode.C
 * Recode beta orbits to 2-adic values.
 * This is the "compressor" function, turned on it's side, across all beta's.
 *
 * Jan 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

// This is described in the text as the "cpr" function, or the
// "compressor function". It just recodes a beta bit expansion
// as a binary bit expansion. Figure 8 or thereabouts, early
// in the diary.
double cpr(double beta, double x)
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
	return cpr(y, x);
}

DECL_MAKE_HEIGHT(coding)
