/*
 * euler-q.C
 *
 * FUNCTION:
 * Display euler q-series aka dedekind eta, 
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


static double euler_prod (double re_q, double im_q)
{
	double rep, imp;
	euler_prod_c (re_q, im_q, &rep, &imp);
	return sqrt (rep*rep+imp*imp);
}


static double dedekind_eta (double re_q, double im_q)
{
	double rep, imp;
	dedekind_eta_c (re_q, im_q, &rep, &imp);
	return sqrt (rep*rep+imp*imp);
}

static double discriminant (double re_q, double im_q)
{
	double rep, imp;
	discriminant_c (re_q, im_q, &rep, &imp);

	// return sqrt (rep*rep+imp*imp);
	return rep;
	// return imp;
}

/*-------------------------------------------------------------------*/
/* This routine fills in the interior of the the convergent area of the 
 * Euler q-series (dedekind eta function) in a simple way 
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
   
   im_position = im_start;
   for (i=0; i<sizey; i++) 
	{
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) 
		{

			double phi = euler_prod (re_position, im_position);
			// double phi = dedekind_eta (re_position, im_position);
			// double phi = discriminant (re_position, im_position);
         glob [i*sizex +j] = phi;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
