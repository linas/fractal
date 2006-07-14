/*
 * divisor.C
 *
 * FUNCTION:
 * display distribution of divisor
 *
 * HISTORY:
 * quick hack -- Linas Vepstas October 1989
 * modernize -- Linas Vepstas March 1996
 * more stuff -- January 2000
 * more stuff -- October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"
#include "harmonic.h"
#include "moebius.h"

/*-------------------------------------------------------------------*/
/* Bifurcation diagram callback, does one row at a time */

static int ex=1;
static long double sum=0.0L;

static void 
divisor_diagram 
(float *array, 
	int array_size, 
	double x_center,
	double x_width,
	double K, 
	int itermax,
	double omega)
{
	int j, cnt=0;

	/* clear out the row */   
	for (j=0; j<array_size; j++)
	{
		array[j] = 0.0;
	}

	for (j=0; j<itermax; j++)
	{
		/* compute divisor summatory function */
		int d = divisor(ex);
		sum += d;

		/* subtract leading term of the function */
		long double del = sum - ex*(logl(ex) + 2.0*M_GAMMA -1.0L);

		/* scale as an estimate of vert size */
		long double scale = 2.0L * powl (ex, 7.0L/22.0L);
		del /= scale;
		del = 0.5L * (del+1.0L);
		ex++;

		/* convert coordinate to integer, and historgram it. */
		double en = array_size * del;
		int n = (int) en;
		if (0 > n) n = 0;
		if (n >= array_size) n = array_size-1;
		array[n] += 1.0;
		cnt ++;
	}
	
	for (j=0; j<array_size; j++)
	{
		array[j] /= cnt;
	}
}


DECL_MAKE_BIFUR(divisor_diagram)

/* --------------------------- END OF LIFE ------------------------- */
