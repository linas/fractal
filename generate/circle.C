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
/* 
 * This routine computes average winding number taken by
 * circle map iterator.
 */
static double rms_winding_number (double omega, double K, int itermax)
{
   double	x=0.0, sq=0.0;
   int		iter,j;
	int cnt=0;
	double start=0.0, end=0.0;
   
	for (j=0; j<itermax/SAMP; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		t -= 0.5;
		x = t;
		start = x;
  
   	/* OK, now start iterating the circle map */
   	for (iter=0; iter < SAMP; iter++) {
      	x += omega - K * sin (2.0 * M_PI * x);
			sq += (x-t)*(x-t);
			t = x;
			cnt ++;
   	}
		end += x;
	}
	
   // x = sqrt (sq) / (end-start);
   x = sqrt (sq) / ((double) cnt);
	// if (K != 0.0) x /= K;
	return x;
}

/*-------------------------------------------------------------------*/
/*
 * Compute the poincare recurrance time for the circle map
 */

#define EPSILON  	0.003
#define SETTLE_TIME 	90
#define RSAMP 200

double 
circle_poincare_recurrance_time (double omega, double K, int itermax)

{
   double	x, y;
   double	xpoint;
   int		j, iter;
   long		num_recurs, time_recur=0;
   
  	num_recurs = 0;
	for (j=0; j<itermax/RSAMP; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		t -= 0.5;
		x = t;
  
   	/* First, we give a spin for 500 cycles, giving the non-chaotic 
    	* parts a chance to phase-lock */
   	for (iter=0; iter<SETTLE_TIME; iter++) 
		{
     		x += omega - K * sin (2.0 * M_PI * x);
   	}
	
   	/* OK, now, we begin to measure the average amount of time to recur */
   	/* (note that we don't have todo += with iter, since its already a running sum). */
   	xpoint = x;
		long ptime = 0;
   	for (iter=0; iter < RSAMP; iter++)
		{
      	x += omega - K * sin (2.0 * M_PI * x);
      	y = fabs (x-xpoint);
      	y -= floor (y);
      	if (y < EPSILON) 
			{
         	num_recurs ++;
         	ptime = iter;
      	}
   	}
		time_recur += ptime;
	}

   /* x is the (normalized) number of cycles to reach recurrance */
   x = (double) time_recur / ((double)num_recurs);

   return x;
}

static double circle_map (double omega, double K, int itermax, double param)
{
	// return winding_number (omega,K, itermax);
	// return noisy_winding_number (omega,K, itermax, param);
	return rms_winding_number (omega,K, itermax);
	// return circle_poincare_recurrance_time (omega,K, itermax);
}

DECL_MAKE_HISTO (circle_map);

/* --------------------------- END OF LIFE ------------------------- */

