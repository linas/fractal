/*
 * euler-q.C
 *
 * FUNCTION:
 * display euler q-series
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

/*-------------------------------------------------------------------*/
/* this routine fills in the interior of the mandelbrot set using */
/* the classic algorithm */

void mandelbrot_stop (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   int		itermax)
{
   int		i,j, globlen;
   double	re_start, im_start, delta;
   double	re_position, im_position;
   double	re, im, tmp;
   int		loop;
   
   delta = width / (double) sizex;
   re_start = re_center - width / 2.0;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;
   
   im_position = im_start;
   for (i=0; i<sizey; i++) {
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) {
         re = re_position;
         im = im_position;
         for (loop=1; loop <itermax; loop++) {
            tmp = re*re - im*im + re_position;
            im = 2.0*re*im + im_position;
            re = tmp;
            if ((re*re + im*im) > 154.0) break;
         }    
         /* glob [i*sizex +j] = (2.0+ re)/2.5;  */
         /* glob [i*sizex +j] = re; */
         glob [i*sizex +j] = im;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
