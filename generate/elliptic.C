/*
 * elliptic.C
 *
 * FUNCTION:
 * display Jacobian elliptic functions, as function of nome,
 * on the unit disk.
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
static double vee = 3.0;

static void elliptic_sn_c (double re_q, double im_q, double *prep, double *pimp)
{
	double tmp;
	int i;
	*prep = 0.0;
	*pimp = 0.0;

	double rep = 0.0;
	double imp = 0.0;


	double qpr = 1.0;
	double qpi = 0.0;

	double qpmod = re_q*re_q+im_q*im_q;
	if (1.0 <= qpmod) return;

   // compute square root of q
	qpmod = sqrt (qpmod);
	double ph = 0.5 * atan2 (im_q, re_q);
   qpr = qpmod * cos (ph);
   qpi = qpmod * sin (ph);

	for (i=0; i<max_terms; i++)
	{
		double term = 2*i+1;
		term = sin (term *vee);

		/* compute (q^(n+1/2))^2 */
		double resq = qpr*qpr - qpi*qpi;
		double imsq = 2.0*qpr*qpi; 
		
		/* compute 1-q^(2n+1) */
		resq = 1.0-resq;

		/* compute 1/(1-q^(2n+1)) */
		qpmod = resq*resq + imsq*imsq;
		resq = resq / qpmod;
		imsq = - imsq / qpmod;

		/* compute full fraction */
		double fre = resq*qpr - imsq * qpi;
		double fim = resq*qpi + imsq * qpr;

		/* accumulate */
		rep += fre *term;
		imp += fim *term;

		/* compute q^k */
		tmp = qpr*re_q - qpi * im_q;
		qpi = qpr*im_q + qpi * re_q;
		qpr = tmp;

		qpmod = qpr*qpr + qpi*qpi;
		if (qpmod < 1.0e-30) break;
	}
	if (max_terms-1 < i)
	{
		// printf ("not converged re=%g im=%g modulus=%g\n", re_q, im_q, qpmod);
	}

	/* multiply by (1-|q|)^2 */
	// tmp = 1.0 - sqrt (re_q*re_q + im_q*im_q);
	// tmp *= tmp;

	// rep *= tmp;
	// imp *= tmp;

	*prep = rep;
	*pimp = imp;
}

static double elliptic_sn (double re_q, double im_q)
{
	double rep, imp;
	elliptic_sn_c (re_q, im_q, &rep, &imp);
	return sqrt (rep*rep+imp*imp);
	// return rep;
}

/*-------------------------------------------------------------------*/
/* This routine fills in the interior of the the convergent area of the 
 * Euler totient in a simple way 
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

	max_terms = itermax;
   
   im_position = im_start;
   for (i=0; i<sizey; i++) 
	{
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) 
		{

			double phi = elliptic_sn (re_position, im_position);
         glob [i*sizex +j] = phi;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
