/*
 * brat-gap-hair.C
 *
 * FUNCTION:
 * Explore Hausdorf measure of mandelbrot set.
 * And other stuff.
 *
 * HISTORY:
 * quick hack -- Linas Vepstas October 1989
 * modernize -- Linas Vepstas March 1996
 * more stuff -- January 2000
 * more stuff -- October 2004
 */

#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "brat.h"
#include "Farey.h"
#include "FareyTree.h"


/*-------------------------------------------------------------------*/
/* The gap-tongue tries to draw the basic gaps
 * of the continued fraction. Te resulting plots look like wavey hair
 * or wavey seaweed. 
 */


void 
MakeHisto (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax,
   double 	renorm)
{
   int		i,j;

   int globlen = sizex*sizey;
   for (i=0; i<globlen; i++) {
      glob [i] = 0.0;
   }

	int d,n;  // denom, numerator
	ContinuedFraction f;

	d = itermax;
	
	for (n=1; n<d; n++)
	{
		int nn = n;
		int dd = d;

// #define DO_TONG_RAND
#ifdef DO_TONG_RAND
		nn = rand() >> 10;
		dd = rand() >> 10;
#endif
		nn %= dd;
		if (0 == nn) continue;
		if (0 == dd) continue;

		int gcf = gcf32 (nn,dd);

		nn /= gcf;
		dd /= gcf;
			
		f.SetRatio (nn,dd);

		double x = (double)nn/(double) dd;

		for (j=0; j<sizey; j++)
		{
			double w = (((double) (sizey-j))-0.5)/((double) sizey);

			// double gap = f.ToXEven(w);
			// double gap = f.ToXPlus(w);
#define REMOVE_LEADING_TERMS
#ifdef REMOVE_LEADING_TERMS
			double gap = f.ToXEven(w);
			gap -= x;
			gap *= (double) dd;
			gap *= (double) dd;
			gap -= w;
			gap += w*w;
			gap -= 0.5*w*w*w;
			gap += x;
#endif
#ifdef WHATEVER
			double gap = f.ToXOdd(w);
			gap -= (1.0-x);
			gap *= (double) dd;
			gap *= (double) dd;
			gap += w;
			gap += w*w;
			gap -= 0.5*w*w*w;
			gap += (1.0-x);
#endif

	// printf ("duude x=%d/%d = %g w=%g gap=%g\n", nn, dd, x, w, gap);

			i = (int) (gap * (double) sizex);
			if (0>i) continue;
			if (i>=sizex) continue;

			glob [j*sizex +i] ++;
		}
	}

   /* renormalize */
	double r = ((double) sizex) / ((double) itermax);
   for (i=0; i<sizex*sizey; i++) 
	{
		glob [i] *= r;
   }
}

/* --------------------------- END OF LIFE ------------------------- */
