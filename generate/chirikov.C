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

#define IMULT 1000
// #define SAMP 1000
#define SAMP 1
	double norm = ((double) sizex*sizey);
	norm /= ((double) SAMP);
	norm /= ((double) itermax);
	norm /= ((double) IMULT);
	norm /= width*height;
	for (k=0; k<SAMP; k++)
	{

  		/* OK, now start iterating the circle map */
		int iter, itermult;
		double pos, posprev, moment;
		pos = rand();
		pos /= RAND_MAX;
		moment = rand();
		moment /= RAND_MAX;
		posprev = 0.0;
  		for (itermult=0; itermult < IMULT; itermult++)
  		for (iter=0; iter < itermax; iter++)
		{
			posprev = pos;
     		pos += moment + K * sin (2.0 * M_PI * pos);
			pos -= floor(pos);

			moment = pos - posprev;
			moment -= floor(moment);

			/* convert to piposel coords */
			double sx = 0.5 + (pos - re_center) / width;
			double sy = 0.5 + (moment - im_center) / height;
			i = (sizex-1) * sx;
			j = (sizey-1) * sy;
			if ((0 <= i) && (i < sizex) &&
			    (0 <= j) && (j < sizey))
			{
				j = sizey-1 - j;
				glob [i + j*sizex] += norm;
			}
  		}
	}
}

/* --------------------------- END OF LIFE ------------------------- */
