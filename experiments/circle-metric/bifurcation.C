/*
 * bifurcation.C
 *
 * FUNCTION:
 * Bifurcation diagram for the standard circle map.
 *
 * HISTORY:
 * quick hack -- Linas Vepstas November 2023
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/*-------------------------------------------------------------------*/
/* Bifurcation diagram callback, does one row at a time */

static void
bifurcation_diagram
(float *array,
	int array_size,
	double x_center,
	double x_width,
	double K,
	int itermax,
	double omega)
{
	long cnt = 0;

	/* clear out the row */
	for (int j=0; j<array_size; j++)
	{
		array[j] = 0.0;
	}

#define ITER_DEPTH 1500  // Number of iteration steps.
	for (int j=0; j<itermax/ITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		/* Settle */
#define SETTLE_TIME 90
		for (int iter=0; iter < SETTLE_TIME; iter++) {
			x += omega - K * sin (2.0 * M_PI * x);
		}

		/* OK, Now iterate the circle map */
		for (int iter=0; iter < ITER_DEPTH; iter++) {
			x += omega - K * sin (2.0 * M_PI * x);

			/* convert coordinate to integer, and historgram it. */
			double y = x - floor(x);
			double en = array_size * y;
			int n = (int) floor(en + 0.5);
			if (0 > n) n = 0;
			if (n >= array_size) n = array_size-1;
			array[n] += 1.0;
			cnt ++;
		}
	}

	// Renormalize
	for (int j=0; j<array_size; j++)
	{
		array[j] *= array_size / ((double) cnt);
	}
}

DECL_MAKE_BIFUR(bifurcation_diagram)

/* --------------------------- END OF LIFE ------------------------- */
