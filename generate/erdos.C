/*
 * erdos.C
 *
 * FUNCTION:
 * Display lambert series involving divisors
 * More generally, 
 * show Weierstrass elliptic function g_2 and g_3 invariants
 * Also show the modular discriminant.
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
#include "modular.h"

static long double erdos_series (long double re_q, long double im_q)
{
	long double rep, imp;
	erdos_series_c (re_q, im_q, 0, &rep, &imp);
	// return sqrt (rep*rep+imp*imp);
	// return imp;
	long double phase = atan2 (imp, rep);
	phase += M_PI;
	phase /= 2.0*M_PI;
	return phase;
}


static long double discriminant (long double re_q, long double im_q)
{
	long double rep, imp;
	disc_c (re_q, im_q, &rep, &imp);

	// return sqrt (rep*rep+imp*imp);
	return rep;
	// return imp;
	// long double phase = atan2 (imp, rep);
	// phase += M_PI;
	// phase /= 2.0*M_PI;
	// return phase;
}

/* Weierstrass elliptic invarient g_2, where q is the nome */
static long double gee_2 (long double re_q, long double im_q)
{
	long double rep, imp;
	gee_2_c (re_q, im_q, &rep, &imp);

	// return sqrt (rep*rep+imp*imp);
	// return rep;
	// return imp;
	long double phase = atan2 (imp, rep);
	phase += M_PI;
	phase /= 2.0*M_PI;
	return phase;
}

static long double gee_3 (long double re_q, long double im_q)
{
	long double rep, imp;
	gee_3_c (re_q, im_q, &rep, &imp);

	// return sqrt (rep*rep+imp*imp);
	return rep;
	// return imp;
	// long double phase = atan2 (imp, rep);
	// phase += M_PI;
	// phase /= 2.0*M_PI;
	// return phase;
}

/*-------------------------------------------------------------------*/
/* This routine fills in the interior of the the convergent area of the 
 * Euler erdos in a simple way 
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
   int		i,j, globlen;
   double	re_start, im_start, delta;
   double	re_position, im_position;
   
   delta = width / (double) sizex;
   re_start = re_center - width / 2.0;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;

	modular_max_terms = itermax;
   
   im_position = im_start;
   for (i=0; i<sizey; i++) 
	{
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) 
		{

			// double phi = erdos_series (re_position, im_position);
			double phi = gee_2 (re_position, im_position);
			// double phi = gee_3 (re_position, im_position);
			// double phi = discriminant (re_position, im_position);
         glob [i*sizex +j] = phi;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
