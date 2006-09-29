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
#include "gcf32.h"

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
   double	re_position, im_position;
   
   delta = width / (double) sizex;
   re_start = re_center - width / 2.0;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = -1.0;

	int p,q;
	for(q=3; q<2000; q+=2)
	{
		for (p=1; p<q; p+=2)
		{
			/* Find fraction p/q reduced to lowest terms */
			int m = gcf32 (p,q);
			if (1<m) continue;

			/* prime the pump */
			int n;
			for (n=0; n<10; n++)
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
		}
	}
         glob [i*sizex +j] = phi;
   }
}

/* --------------------------- END OF LIFE ------------------------- */
