/*
 * circle.C
 *
 * FUNCTION:
 * Circle map -- do it all over again.
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

static int max_terms;

static double winding_number (double omega, double K, int itermax)
{
	
}

static double circle_map (double omega, double K, int itermax)
{
	return winding_number (omega,K, itermax);
}

DECL_MAKE_HISTO (circle_map);

/* --------------------------- END OF LIFE ------------------------- */

/*-------------------------------------------------------------------*/
/* 
 * This routine computes average winding number taken by
 * circle map iterator.
 */

float circle_winding_number (CircleData *dat,
                       double omega,
                       double K)

{
   double	x=0.0;
   int		iter;
   
  
   /* OK, now start iterating the circle map */
   for (iter=0; iter < dat->itermax; iter++) {
      x += omega - K * sin (2.0 * M_PI * x);
   }

   /* x is the normalized winding number */
   x /= (double) (dat->itermax);

/* ============
   if (x<0.0) x = 0.0;
   if (x>1.0) x = 1.0;
   if ((x<0.0) | (x>1.0)) {
      printf ("out of bounds omega=%f and K=%f W=%f \n", omega, K, x);
   }
   printf ("for omega=%f and K=%f W=%f \n", omega, K, x);
   ==============*/

   return ((float) x);
}
   
