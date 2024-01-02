/*
 * point.C
 *
 * Dec 2017 Jan 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/*-------------------------------------------------------------------*/
/*
 */

static void bifurcation_diagram (float *array,
                                 int array_size,
                                 double x_center,
                                 double x_width,
                                 double beta,
                                 int itermax,
                                 double omega)
{
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	// printf("beta=%g\n", beta);

	double midpnt = 0.5*beta;
	double obn = 1.0;
	for (int j=0; j<500; j++)
	{
		double x = midpnt;

		double en = array_size * (x-floor(x));
		int n = en;
		if (0 > n) n = 0;
		if (n >= array_size) n = array_size-1;
		array[n] += obn;

		if (0.5 < midpnt) midpnt -= 0.5;
		midpnt *= beta;

		obn /= beta;
		if (obn < 1.0e-14) break;
	}

#if 0
	double norm = ((double) array_size) / ((double) cnt);
	for (int j=0; j<array_size; j++)
		array[j] *= norm;
#endif

}

DECL_MAKE_BIFUR(bifurcation_diagram)

#if 0
int main ( int argc, char * argv[])
{
	double om = atof(argv[1]);
	double kb = atof(argv[2]);
}
#endif
