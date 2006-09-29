/*
 * ising-moment.C
 *
 * FUNCTION:
 * Display wave-vector of a periodic point
 *
 * HISTORY:
 * quick hack -- Linas Vepstas Sept 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"
#include "gcf.h"

/*-------------------------------------------------------------------*/
/* This routine fills in the interior of the Baker's map square */

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
   int		i,j, globlen;
   double	re_start, im_start, delta;
   
   delta = width / (double) sizex;
   re_start = re_center - width / 2.0;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = -1.0;

	int p,q;
	for(q=3; q<itermax; q+=2)
	{
		for (p=1; p<q; p+=2)
		{
			/* Find fraction p/q reduced to lowest terms */
			int m = gcf32 (p,q);
			if (1<m) continue;

			// double frac = ((double) p)/((double) q);
			double frac = 1.0/((double) q);

			/* iterate past any decaying components */
			int n;
			for (n=0; n<15; n++)
			{
				p *= 2;
				if (p < q)
				{
					/* bit is zero */
				}
				else
				{
					/* bit is one */
					p -= q;
				}
			}
			if (0 == p) continue;

#if 1
			/* measure the period */
			int pstart = p;
			frac = 0;
			for (n=0; n<q; n++)
			{
				p *= 2;
				if (p < q)
				{
					/* bit is zero */
				}
				else
				{
					/* bit is one */
					p -= q;
				}
				frac ++;
				if (p == pstart) break;
			}
#endif

#define GRID 9
			/* prime the pump -- fill bits left to right,
			 * with j on the left and i on the right. */
			i=j=0;
			for (n=0; n<GRID; n++)
			{
				p *= 2;
				if (p < q)
				{
					/* bit is zero */
				}
				else
				{
					/* bit is one */
					j += 1<<n;
					p -= q;
				}
			}
			for (n=0; n<GRID; n++)
			{
				p *= 2;
				if (p < q)
				{
					/* bit is zero */
				}
				else
				{
					/* bit is one */
					i += 1<<(GRID-1-n);
					p -= q;
				}
			}

			/* now fill the array */
			for (n=0; n<60; n++)
			{
				if (0.0 > glob [i*sizex +j])
				{
					glob [i*sizex +j] = frac;
				}

				/* shift the bits over to the left */
				j /= 2;
				if (i & (1<<(GRID-1)))
				{
					j += 1<<(GRID-1);
					i -= 1<<(GRID-1);
				}
				i *= 2;

				/* and add a bit on the far right */
				p *= 2;
				if (p < q)
				{
					/* bit is zero */
				}
				else
				{
					/* bit is one */
					i += 1;
					p -= q;
				}
			}
		}
	}
}

/* --------------------------- END OF LIFE ------------------------- */
