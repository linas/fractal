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

/* Return euler-product form of the q-series (dedekind eta) */

static void euler_prod_c (double re_q, double im_q, double *prep, double *pimp)
{
	int i;
	*prep = 0.0;
	*pimp = 0.0;

	double rep = 1.0;
	double imp = 0.0;

	double qpr = re_q;
	double qpi = im_q;

	double qpmod = qpr*qpr+qpi*qpi;
	if (1.0 <= qpmod) return;

	for (i=0; i<100100; i++)
	{
		double tmp;

		/* compute prod (1-q^k) rember to use 1-real and -imag .. */
		tmp = (1.0-qpr)*rep + qpi*imp;
		imp = (1.0-qpr)*imp - qpi*rep;
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

	*prep = rep;
	*pimp = imp;
}

static double euler_prod (double re_q, double im_q)
{
	double rep, imp;
	euler_prod_c (re_q, im_q, &rep, &imp);
	return sqrt (rep*rep+imp*imp);
}

/* The dedekind eta multiplies by an additonal factor of q^1/24 */
static void dedekind_eta_c (double re_q, double im_q, double *pre, double *pim)
{
	double rep, imp;
	euler_prod_c (re_q, im_q, &rep, &imp);

	double phase = atan2 (im_q, re_q);
	phase /= 24.0;
	double mod = sqrt (re_q*re_q + im_q*im_q);
	mod = pow (mod, 1.0/24.0);

	double reqt = mod * cos (phase);
	double imqt = mod * sin (phase);

	double tmp = reqt*rep - imqt*imp;
	imp = reqt*imp + imqt*rep;
	rep = tmp;

	*pre = rep;
	*pim = imp;
}

static double dedekind_eta (double re_q, double im_q)
{
	double rep, imp;
	dedekind_eta_c (re_q, im_q, &rep, &imp);
	return sqrt (rep*rep+imp*imp);
}

/* modular discriminat = eta to 24 */
static double discriminant (double re_q, double im_q)
{
	double rep, imp;
	euler_prod_c (re_q, im_q, &rep, &imp);

	/* take euler product to 24'th power */
	double phase = atan2 (imp, rep);
	phase *= 24.0;
	double mod = sqrt (rep*rep + imp*imp);
	mod = pow (mod, 24.0);

	double red = mod * cos (phase);
	double imd = mod * sin (phase);

	/* multply by q one more time .. */
	double tmp = red * re_q - imd * im_q;
	imd = red * im_q + imd * re_q;
	red = tmp;

	/* multiply by 2pi to the 12'th */
	mod = pow (2.0*M_PI, 12.0);
	red *= mod;
	imd *= mod;


	// return sqrt (red*red+imd*imd);
	return red;
	// return imd;
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

			// double phi = euler_prod (re_position, im_position);
			// double phi = dedekind_eta (re_position, im_position);
			double phi = discriminant (re_position, im_position);
         glob [i*sizex +j] = phi;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
