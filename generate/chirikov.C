/*
 * chirikov.C
 *
 * FUNCTION:
 * Chirikov-Taylor Standard Map
 *
 * HISTORY:
 * quick hack -- Linas Vepstas November 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/*-------------------------------------------------------------------*/
/*
 * This routine computes a scatterplot of position vs momentum for
 * the standard map.
 */

void MakeHisto (
   float    *glob,
   int      sizex,
   int      sizey,
   double   re_center,
   double   im_center,
   double   width,
   double   height,
   int      itermax,
   double   param)
{
	int      i,j, k, globlen;

	globlen = sizex*sizey;
	for (i=0; i<globlen; i++) glob [i] = 0.0;

	double K = param/ (2.0*M_PI);

#define SAMP 1000
	double norm = ((double) sizex*sizey) / ((double) SAMP*itermax);
	for (k=0; k<SAMP; k++)
	{

  		/* OK, now start iterating the circle map */
		int iter;
		double x, xprev, moment;
		x = rand();
		x /= RAND_MAX;
		moment = rand();
		moment /= RAND_MAX;
		xprev = 0.0;
  		for (iter=0; iter < itermax; iter++)
		{
			xprev = x;
     		x += moment + K * sin (2.0 * M_PI * x);
			x -= floor(x);

			moment = x - xprev;
			moment -= floor(moment);

			/* convert to pixel coords */
			i = (sizex-1) * x;
			j = (sizey-1) * moment;
			j = sizey-1 - j;

			glob [i + j*sizex] += norm;
  		}
	}
}

/* --------------------------- END OF LIFE ------------------------- */
