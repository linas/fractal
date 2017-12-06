/*
 * bernie.C
 *
 * Altered simplified algorithmic Bernoulli map
 * Dec 2017
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bitops.h"
#include "brat.h"

/*-------------------------------------------------------------------*/
/*
 */

static double xprod(double x, double y, int itermax, double param)
{
	// return mult_xor(x, y);
	return x*y;
}

DECL_MAKE_HEIGHT (xprod);

#if 0
int main ( int argc, char * argv[])
{
	double om = atof(argv[1]);
	double kb = atof(argv[2]);
}
#endif
