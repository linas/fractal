/*
 * mobius.C
 *
 * FUNCTION:
 * display mobius mu function
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
#include "moebius.h"

int thue_morse(int n)
{
	if (0 == n) return 0;
	if (1 == n) return 1;
	if (0 == n%2) return thue_morse (n/2);
	return (1-thue_morse ((n-1)/2));
}

double randoid(int n)
{
#define NVAL 65186
	static double array[NVAL];
	static int inited=0;
	if (!inited)
	{
		int i;
		inited = 1;
		srand (99);
		for (i=0; i<NVAL; i++)
		{
			// array[i] = rand() & 0x1;
			array[i] = ((double) rand()) / ((double)RAND_MAX);
		}
	}

	if (NVAL<= n) return 0;
	return array[n];
}

static int max_terms;

static void mobius_series_c (double re_q, double im_q, double *prep, double *pimp)
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

	for (i=0; i<max_terms; i++)
	{

		// double t = moebius_mu (i+1);
		// double t = mertens_m (i+1);
		// double t = liouville_omega (i+1);
		// double t = liouville_lambda (i+1);
		// double t = mangoldt_lambda (i+1);
		// double t = thue_morse (i+1);
		// int tm = thue_morse (i+1);
		// double t = 1.0;
		// if (1 == tm) t = -1.0;

		// double t = moebius_mu (i+1);
		double t = randoid (i+1);
		t *= (i+1);
		t *= (i+1);
#if 0
		t *= (i+1);
#endif

		rep += qpr *t;
		imp += qpi *t;

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

	*prep = rep;
	*pimp = imp;
}

static void curvature_c (double re_q, double im_q, double *prep, double *pimp)
{
	double tmp;
	int i;
	*prep = 0.0;
	*pimp = 0.0;

	double rep = 0.0;
	double imp = 0.0;

	double qpr = 1.0;
	double qpi = 0.0;

	double qprm1 = 0.0;
	double qpim1 = 0.0;

	double qprm2 = 0.0;
	double qpim2 = 0.0;

	double qpmod = re_q*re_q+im_q*im_q;
	if (1.0 <= qpmod) return;

	for (i=0; i<max_terms; i++)
	{

		// double t = moebius_mu (i+1);
		// double t = mertens_m (i+1);
		// double t = liouville_omega (i+1);
		// double t = liouville_lambda (i+1);
		// double t = mangoldt_lambda (i+1);
		// double t = thue_morse (i+1);
		// int tm = thue_morse (i+1);
		// double t = 1.0;
		// if (1 == tm) t = -1.0;

		// double t = moebius_mu (i+1);
		double t = randoid (i+1);
#if 0
		t *= (i+1);
		t *= (i+1);
		t *= (i+1);
#endif

		double eye = i;
		double rehp = eye*qprm1*t;
		double imhp = eye*qpim1*t;
		double rehpp = eye*(eye-1.0)*qprm2*t;
		double imhpp = eye*(eye-1.0)*qpim2*t;

		double norm = pow (rehp*rehp+imhp*imhp,  1.5);
		double  equipot_term = - rehpp*(rehp*rehp - imhp*imhp) - 2.0*rehp*imhp*imhpp;
		equipot_term /= norm;

		double ray_term = imhpp*(rehp*rehp - imhp*imhp) - 2.0*rehp*imhp*rehpp;
		ray_term /= norm;

		rep += equipot_term;
		imp += ray_term;

		/* save lower derives */
		qprm2 = qprm1;
		qpim2 = qpim1;

		qprm1 = qpr;
		qpim1 = qpi;

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

	*prep = rep;
	*pimp = imp;
}


static double mobius_series (double re_q, double im_q)
{
	double rep, imp;

	// mobius_series_c (re_q, im_q, &rep, &imp);
	curvature_c (re_q, im_q, &rep, &imp);
	// return sqrt (rep*rep+imp*imp);
	return rep;
	// return imp;
	// return (atan2 (imp,rep)+M_PI)/(2.0*M_PI);
}

/*-------------------------------------------------------------------*/
/* This routine fills in the interior of the the convergent area of the 
 * Euler mobius in a simple way 
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

			double phi = mobius_series (re_position, im_position);
         glob [i*sizex +j] = phi;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
