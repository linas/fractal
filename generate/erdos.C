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
#include "moebius.h"


static int max_terms;

/*
 * Compute the divisor arithmetic function
 */

int divisor (int n)
{
	int acc = 0;
	int d;

	for (d=1; d<=n; d++)
	{
		if (n%d) continue;
		acc ++;
	}

	return acc;
}

/* An erdos-borwein-like series sum_n x^n/(1-x^n) 
 */
static void erdos_series_c (double re_q, double im_q, int sigma, double *prep, double *pimp)
{
	double tmp;
	int i;
	*prep = 0.0;
	*pimp = 0.0;

	double rep = 0.0;
	double imp = 0.0;

	double qpr = 1.0;
	double qpi = 0.0;

	qpr = re_q;
	qpi = im_q;

	double qpmod = re_q*re_q+im_q*im_q;
	if (1.0 <= qpmod) return;

	long double tn = 0.5;
	for (i=0; i<max_terms; i++)
	{
		double dr, di;

#define LAMBERT_SUM
#ifdef LAMBERT_SUM
		/* compute 1/(1-q^n) OK */
		tmp = (1.0-qpr)*(1.0-qpr) + qpi*qpi;
		tmp = 1.0/tmp;
		dr = (1.0-qpr) * tmp;
		di = qpi * tmp;

		/* compute q^n/(1-q^n) */
		tmp = dr*qpr - di*qpi;
		di = dr*qpi + di*qpr;
		dr = tmp;

		/* Multiply by power representing sigma_k 
		 * sigma = 3 for g_2 and 5 for g_3 
		 */
		int j;
		for (j=0; j<sigma; j++) {
			dr *= i+1;
			di *= i+1;
		}

#endif

#ifdef DIVISOR_SUM
		tmp = divisor (i+1);
		dr = qpr*tmp;
		di = qpi*tmp;
#endif

#if 0
		tmp = 1.0L/(1.0L-tn);
		tmp /= i+1;
		dr = qpr*tmp;
		di = qpi*tmp;
#endif

		rep += dr;
		imp += di;

		/* compute q^k */
		tmp = qpr*re_q - qpi * im_q;
		qpi = qpr*im_q + qpi * re_q;
		qpr = tmp;

		qpmod = qpr*qpr + qpi*qpi;
		if (qpmod < 1.0e-30) break;

		tn *= 0.5L;
	}
	if (max_terms-1 < i)
	{
		// printf ("not converged re=%g im=%g modulus=%g\n", re_q, im_q, qpmod);
	}

	*prep = rep;
	*pimp = imp;
}

static double erdos_series (double re_q, double im_q)
{
	double rep, imp;
	erdos_series_c (re_q, im_q, 0, &rep, &imp);
	// return sqrt (rep*rep+imp*imp);
	// return imp;
	double phase = atan2 (imp, rep);
	phase += M_PI;
	phase /= 2.0*M_PI;
	return phase;
}

/* Weierstrass elliptic invarient g_2, where q is the nome */
static void gee_2_c (double re_q, double im_q, double *pre, double *pim)
{
	double rep, imp;
	double sqre, sqim;

	sqre = re_q*re_q - im_q *im_q;
	sqim = 2.0*re_q * im_q;
	
	erdos_series_c (sqre, sqim, 3, &rep, &imp);
	rep *= 240.0;
	imp *= 240.0;
	rep +=1.0;
	rep *= 4.0 *M_PI*M_PI*M_PI*M_PI / 3.0;
	imp *= 4.0 *M_PI*M_PI*M_PI*M_PI / 3.0;

	*pre = rep;
	*pim = imp;
}

static void gee_3_c (double re_q, double im_q, double *pre, double *pim)
{
	double rep, imp;
	double sqre, sqim;

	sqre = re_q*re_q - im_q *im_q;
	sqim = 2.0*re_q * im_q;
	
	erdos_series_c (sqre, sqim, 5, &rep, &imp);
	rep *= -504.0;
	imp *= -504.0;
	rep +=1.0;
	rep *= 8.0 *M_PI*M_PI*M_PI*M_PI *M_PI*M_PI/ 27.0;
	imp *= 8.0 *M_PI*M_PI*M_PI*M_PI *M_PI*M_PI/ 27.0;

	*pre = rep;
	*pim = imp;
}

/* the modular discriminant */
static void disc_c (double re_q, double im_q, double *pre, double *pim)
{
	double g2re, g2im;
	double g3re, g3im;
	double tmp;

	gee_2_c (re_q, im_q, &g2re, &g2im);
	gee_3_c (re_q, im_q, &g3re, &g3im);
	
	double g3sqre, g3sqim;
	g3sqre = g3re*g3re - g3im*g3im;
	g3sqim = 2.0 * g3re*g3im;
	
	double g2cure, g2cuim;
	g2cure = g2re*g2re - g2im*g2im;
	g2cuim = 2.0 * g2re*g2im;
	tmp = g2cure*g2re - g2cuim*g2im;
	g2cuim = g2cure*g2im + g2cuim* g2re;
	g2cure = tmp;

	double dre, dim;
	dre = g2cure - 27.0*g3sqre;
	dim = g2cuim - 27.0*g3sqim;

	*pre = dre;
	*pim = dim;
}

static double discriminant (double re_q, double im_q)
{
	double rep, imp;
	disc_c (re_q, im_q, &rep, &imp);

	// return sqrt (rep*rep+imp*imp);
	return rep;
	// return imp;
	// double phase = atan2 (imp, rep);
	// phase += M_PI;
	// phase /= 2.0*M_PI;
	// return phase;
}

/* Weierstrass elliptic invarient g_2, where q is the nome */
static double gee_2 (double re_q, double im_q)
{
	double rep, imp;
	gee_2_c (re_q, im_q, &rep, &imp);

	// return sqrt (rep*rep+imp*imp);
	// return rep;
	// return imp;
	double phase = atan2 (imp, rep);
	phase += M_PI;
	phase /= 2.0*M_PI;
	return phase;
}

static double gee_3 (double re_q, double im_q)
{
	double rep, imp;
	gee_3_c (re_q, im_q, &rep, &imp);

	// return sqrt (rep*rep+imp*imp);
	// return rep;
	// return imp;
	double phase = atan2 (imp, rep);
	phase += M_PI;
	phase /= 2.0*M_PI;
	return phase;
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

	max_terms = itermax;
   
   im_position = im_start;
   for (i=0; i<sizey; i++) 
	{
      if (i%10==0) printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) 
		{

			// double phi = erdos_series (re_position, im_position);
			// double phi = gee_2 (re_position, im_position);
			// double phi = gee_3 (re_position, im_position);
			double phi = discriminant (re_position, im_position);
         glob [i*sizex +j] = phi;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
