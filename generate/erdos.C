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

/* sigma arithmetic series, equals divisor for a=0 */
int sigma (int n, int a)
{
	int acc = 0;
	int d;

	for (d=1; d<=n; d++)
	{
		if (n%d) continue;

		int dp = 1;
		int ia;
		for (ia=0; ia<a; ia++) dp *= d;
		acc += dp;
	}

	return acc;
}

/* An erdos-borwein-like series sum_n x^n/(1-x^n) 
 */
static void erdos_series_c (long double re_q, 
                            long double im_q, 
                            int sa, long double *prep, long double *pimp)
{
	long double tmp;
	int i;
	*prep = 0.0L;
	*pimp = 0.0L;

	long double rep = 0.0L;
	long double imp = 0.0L;

	long double qpr = 1.0L;
	long double qpi = 0.0L;

	qpr = re_q;
	qpi = im_q;

	long double qpmod = re_q*re_q+im_q*im_q;
	if (1.0L <= qpmod) return;

	long double tn = 0.5L;
	for (i=0; i<max_terms; i++)
	{
		long double dr, di;

#define LAMBERT_SUM
#ifdef LAMBERT_SUM
		/* compute 1/(1-q^n) OK */
		tmp = (1.0L-qpr)*(1.0L-qpr) + qpi*qpi;
		tmp = 1.0L/tmp;
		dr = (1.0L-qpr) * tmp;
		di = qpi * tmp;

		/* compute q^n/(1-q^n) */
		tmp = dr*qpr - di*qpi;
		di = dr*qpi + di*qpr;
		dr = tmp;

		/* Multiply by power representing sigma_k 
		 * sigma = 3 for g_2 and 5 for g_3 
		 */
		int j;
		for (j=0; j<sa; j++) {
			dr *= (long double) i+1;
			di *= (long double) i+1;
		}

#endif

// #define DIVISOR_SUM
#ifdef DIVISOR_SUM
		// tmp = divisor (i+1);
		tmp = sigma (i+1, sa);
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
		if (qpmod < 1.0e-50) break;

		tn *= 0.5L;
	}
	if (max_terms-1 < i)
	{
		// printf ("not converged re=%g im=%g modulus=%g\n", re_q, im_q, qpmod);
	}

	*prep = rep;
	*pimp = imp;
}

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

/* Weierstrass elliptic invarient g_2, where q is the nome */
static void gee_2_c (long double re_q, long double im_q, long double *pre, long double *pim)
{
	long double rep, imp;
	long double sqre, sqim;

	// sqre = re_q*re_q - im_q *im_q;
	// sqim = 2.0*re_q * im_q;
	sqre = re_q;
	sqim = im_q;
	
	erdos_series_c (sqre, sqim, 3, &rep, &imp);
	rep *= 240.0;
	imp *= 240.0;
	rep +=1.0;
	rep *= 4.0 *M_PI*M_PI*M_PI*M_PI / 3.0;
	imp *= 4.0 *M_PI*M_PI*M_PI*M_PI / 3.0;

	*pre = rep;
	*pim = imp;
}

static void gee_3_c (long double re_q, long double im_q, long double *pre, long double *pim)
{
	long double rep, imp;
	long double sqre, sqim;

	// sqre = re_q*re_q - im_q *im_q;
	// sqim = 2.0*re_q * im_q;
	sqre = re_q;
	sqim = im_q;
	
	erdos_series_c (sqre, sqim, 5, &rep, &imp);
	rep *= -504.0L;
	imp *= -504.0L;
	rep +=1.0L;
	rep *= 8.0L *M_PI*M_PI*M_PI*M_PI *M_PI*M_PI/ 27.0L;
	imp *= 8.0L *M_PI*M_PI*M_PI*M_PI *M_PI*M_PI/ 27.0L;

	*pre = rep;
	*pim = imp;
}

/* the modular discriminant */
static void disc_c (long double re_q, long double im_q, long double *pre, long double *pim)
{
	long double g2re, g2im;
	long double g3re, g3im;
	long double tmp;

	gee_2_c (re_q, im_q, &g2re, &g2im);
	gee_3_c (re_q, im_q, &g3re, &g3im);
	
	long double g3sqre, g3sqim;
	g3sqre = g3re*g3re - g3im*g3im;
	g3sqim = 2.0L * g3re*g3im;
	
	long double g2cure, g2cuim;
	g2cure = g2re*g2re - g2im*g2im;
	g2cuim = 2.0L * g2re*g2im;
	tmp = g2cure*g2re - g2cuim*g2im;
	g2cuim = g2cure*g2im + g2cuim* g2re;
	g2cure = tmp;

	long double dre, dim;
	dre = g2cure - 27.0L*g3sqre;
	dim = g2cuim - 27.0L*g3sqim;

// if(dre < 1.0e-10 *g2cure) printf ("duude bad converge for %g %g == %g\n", re_q, im_q, dre);
// if(fabsl(dre) < 1.0e-10 *fabsl(g2cure)) dre = dim = 1.0;

	*pre = dre;
	*pim = dim;
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

	max_terms = itermax;
   
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
