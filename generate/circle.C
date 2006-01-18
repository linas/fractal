/*
 * circle.C
 *
 * FUNCTION:
 * Circle map -- all new, all improved -- do it all over again.
 *
 * HISTORY:
 * quick hack -- Linas Vepstas October 1989
 * modernize -- Linas Vepstas March 1996
 * more stuff -- January 2000
 * more stuff -- October 2004
 * more stuff -- January 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/*-------------------------------------------------------------------*/
/* 
 * This routine computes average winding number taken by
 * circle map iterator.
 */
static double winding_number (double omega, double K, int itermax)
{
   double	x=0.0;
   int		iter,j;
	int cnt=0;
	double start=0.0, end=0.0;
   
#define SAMP 100
	for (j=0; j<itermax/SAMP; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		t -= 0.5;
		x = t;
		start += x;
  
   	/* OK, now start iterating the circle map */
   	for (iter=0; iter < SAMP; iter++) {
      	x += omega - K * sin (2.0 * M_PI * x);
			cnt ++;
   	}
		end += x;
	}
	
   x = (end-start) / ((double) cnt);
	return x;
}

/*-------------------------------------------------------------------*/
/* 
 * This routine computes average winding number taken by
 * circle map iterator.
 */
static double rms_winding_number (double omega, double K, int itermax)
{
   double	x=0.0, sq=0.0;
   int		iter,j;
	int cnt=0;
   
#define SAMP 100
	for (j=0; j<itermax/SAMP; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		t -= 0.5;
		x = t;
  
   	/* OK, now start iterating the circle map */
   	for (iter=0; iter < SAMP; iter++) {
      	x += omega - K * sin (2.0 * M_PI * x);
			sq += (x-t-omega)*(x-t-omega);
			t = x;
			cnt ++;
   	}
	}
	
   x = sqrt (sq) / ((double) cnt);
	if (K != 0.0) x /= K;
	return x;
}

static double circle_map (double omega, double K, int itermax)
{
	// return winding_number (omega,K, itermax);
	return rms_winding_number (omega,K, itermax);
}

DECL_MAKE_HISTO (circle_map);

/* --------------------------- END OF LIFE ------------------------- */

