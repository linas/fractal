/*
 * matrix.C
 *
 * Matrix elements visualization
 * Jan 2018
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/*-------------------------------------------------------------------*/
/*
 */
#define NOMAIN 1
#include "psi.c"

static void matrix_diagram (float *array,
                             int array_size,
                             double x_center,
                             double x_width,
                             double K,
                             int itermax,
                             double omega)
{
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;
	find_midpoints(omega);

	int i = K;
	i = 48 - i;
	printf("working i=%d omega=%g\n", i, omega);

	int js = i-1;
	if (js < 0) js = 0;
	for (int j=js; j<array_size; j++)
	{
		array[j] = hess(omega, i, j);
	}
}

DECL_MAKE_BIFUR(matrix_diagram)
