/*
 * borel.C
 *
 * Draw diagram of where thhe midpoint is visiting.
 *
 * Jan 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

#define NOMAIN
#include "visitation.c"

/*-------------------------------------------------------------------*/

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

	visitation_tree(beta);

	int ja = 0;
	double xdelt = 1.0 / ((double) array_size);
	double ex = 0.0;

	for (int i=0; i<DEPTH; i++)
	{
		// Get the x,y coordinates
		if (NEG_ONE == visit[i]) break;
		int j = sorted[i];
		unsigned long v = visit[j];
		double dya = canonical_dyadic(v);
		double midpnt = midp[j];
		midpnt /= 0.5*beta;

		// Copy them into array.
		while (ex < dya && ja < array_size)
		{
			array[ja] = midpnt;
			ja++;
			ex += xdelt;
		}

		if (array_size <= ja) break;
	}

	float last = array[ja-1];
	while (ja < array_size)
	{
		array[ja] = last;
		ja++;
	}
}

DECL_MAKE_BIFUR(bifurcation_diagram)
