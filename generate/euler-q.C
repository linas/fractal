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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/* Return euler-product form of the q-series (dedekind eta) */

static double euler_prod (double re_q, double im_q)
{
	int i;

	double rep = 1.0;
	double imp = 1.0;

	double qpr = re_q;
	double qpi = im_q;

	double qpmod = qpr*qpr+qpi*qpi;
	if (1.0 <= qpmod) return 0.0;

	for (i=0; i<100100; i++)
	{
		double tmp;

		/* compute prod (1-q^k) */
		tmp = (1.0-qpr)*rep - qpi*imp;
		imp = (1.0-qpr)*imp + qpi*rep;
		rep = tmp;

		/* compute q^k */
		tmp = qpr*re_q - qpi * im_q;
		qpi = qpr*im_q + qpi * re_q;
		qpr = tmp;

		qpmod = qpr*qpr + qpi*qpi;
		if (qpmod < 1.0e-30) break;
	}
	if (90100 < i)
	{
		printf ("not converged re=%g im=%g modulus=%g\n", re_q, im_q, qpmod);
	}

	return sqrt (rep*rep+imp*imp);
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
         glob [i*sizex +j] = phi;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
