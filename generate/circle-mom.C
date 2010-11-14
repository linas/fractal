/*
 * circle-mom.C
 *
 * FUNCTION:
 * Circle map momentum.
 *
 * HISTORY:
 * quick hack -- Linas Vepstas October 1989
 * modernize -- Linas Vepstas March 1996
 * more stuff -- January 2000
 * more stuff -- October 2004
 * more stuff -- January 2006
 * momentum graph -- Nov 2010
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
   
#define SAMP 150
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
 * circle map iterator. -- subject to noise
 */
static double noisy_winding_number (double omega, double K, int itermax, double noise)
{
   double	x=0.0;
   int		iter;
	int cnt=0;
   
  	/* OK, now start iterating the circle map */
  	for (iter=0; iter < itermax; iter++) {
     	x += omega - K * sin (2.0 * M_PI * x);
		cnt ++;

		/* white noise, equi-distributed, sharp cutoff */
		double t = rand();
		t /= RAND_MAX;
		t -= 0.5;
		x += noise*t;
  	}
	
   x /= ((double) cnt);
	return x;
}

/*-------------------------------------------------------------------*/

void MakeHisto (
   float    *glob,
   int      sizex,
   int      sizey,
   double   re_center,
   double   im_center,
   double   width,
   double   height,
   int      itermax,
   double   renorm)
{

}

/* --------------------------- END OF LIFE ------------------------- */


