/*
 * circle-mom.C
 *
 * FUNCTION:
 * Circle map momentum.
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
 * the circle map
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

#define SAMP 100
	for (k=0; k<SAMP; k++)
	{
		double omega = 0.333;
		double K = param;

  		/* OK, now start iterating the circle map */
		int iter;
		double x, xprev, moment;
		x = rand();
		x /= RAND_MAX;
		xprev = 0.0;
  		for (iter=0; iter < itermax; iter++)
		{
			xprev = x;
     		x += omega - K * sin (2.0 * M_PI * x);
			x -= floor(x);

			moment = x - xprev;
			moment -= floor(moment);

			/* convert to pixel coords */
			i = (sizex-1) * x;
			j = (sizey-1) * moment;

			glob [i + j*sizex] += 1.0;
  		}
	}
}

/* --------------------------- END OF LIFE ------------------------- */
