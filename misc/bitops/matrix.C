/*
 * matrix.C
 *
 * Visualization of the matrix elements in the Hessenberg basis.
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
#include "psibig.c"

static void matrix_diagram (float *array,
                             int array_size,
                             double x_center,
                             double x_width,
                             double row,
                             int itermax,
                             double K)
{
	static bool init=false;
	if (not init)
	{
		init = true;
		// find_midpoints(K, MAXN);
		big_midpoints(K, 400, midpoints, MAXN);
		sequence_midpoints(K, MAXN);
		printf("working K=%g\n", K);
	}

	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	int i = row;
	i = array_size - i;
	// if (0 == i%20) printf("working i=%d K=%g\n", i, K);

	int js = i-1;
	if (js < 0) js = 0;
	for (int j=js; j<array_size; j++)
	{
		array[j] = hess(K, i, j);
	}
}

DECL_MAKE_BIFUR(matrix_diagram)
