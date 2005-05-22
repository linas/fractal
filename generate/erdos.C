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

void 
poincare_disk_to_plane_coords (double x, double y, 
                               double *px, double *py)
{
	double deno = (1.0-x)*(1.0-x) + y*y;
	deno = 1.0/deno;
	*px = -2.0*y*deno;
	*py = (1.0 - x*x - y*y) * deno;
}

void 
plane_to_poincare_disk_coords (double x, double y, 
                               double *px, double *py)
{
	double deno = x*x + (y+1.0)*(y+1.0);
	deno = 1.0/deno;
	*px = (x*x + y*y - 1.0) * deno;
	*py = -2.0*deno; 
}

void 
plane_to_q_disk_coords (double tau_re, double tau_im, 
                        double *px, double *py)
{
	/* now go back to q-series coords */
	double rq = exp (-tau_im * 2.0 * M_PI);
	*px = rq * cos (tau_re * 2.0 * M_PI);
	*py = rq * sin (tau_re * 2.0 * M_PI);
}

void 
mobius_xform (double a, double b, double c, double d,
              double tau_re, double tau_im, 
              double *px, double *py)
{
	/* now apply mobius */
	double deno = c*tau_re+d;
	deno = deno*deno + c*c*tau_im*tau_im;
	deno = 1.0/deno;

	tau_re = (a*tau_re+b)*(c*tau_re+d) + a*c*tau_re*tau_im;
	*px = tau_re * deno;
	*py = tau_im * deno;
}

long double erdos_series (long double re_q, long double im_q)
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


long double discriminant (long double re_q, long double im_q)
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
long double gee_2 (long double re_q, long double im_q)
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

long double gee_3 (long double re_q, long double im_q)
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

double klein_j (double re_q, double im_q)
{
	double rep, imp;
	klein_j_invariant_c (re_q, im_q, &rep, &imp);

	return sqrt (rep*rep+imp*imp);
	// return rep;
	// return imp;
	// long double phase = atan2 (imp, rep);
	// phase += M_PI;
	// phase /= 2.0*M_PI;
	// return phase;
}

double domain_bounce (double re_t, double im_t, 
                      double (*f)(double,double), int depth)
{
#define DET 5
	depth --;

	double m = re_t*re_t+im_t*im_t;
	double sre = -re_t / m;
	double sim = im_t / m;

	double acc = 0.0;
	int n;
	for (n=-DET; n<=DET; n++)
	{
		double en = n;
		acc += (f) (sre+n, sim);
		if (depth)
			acc += domain_bounce (sre+n, sim, f, depth);
	}

	return acc;
}

double bound (double re_tau, double im_tau)
{
#define BND 0.01
	double m = sqrt (re_tau*re_tau+im_tau*im_tau);
	
	if (fabs(1.0-m) < BND) return 1.0;
	if (m>1.0 && fabs(re_tau-0.5) < BND) return 1.0;
	if (m>1.0 && fabs(re_tau+0.5) < BND) return 1.0;
	return 0.0;

#if 0
	double dist = 0.0;
	dist += BND*fabs (1.0-m);
	dist += BND*fabs (re_tau-0.5);
	dist += BND*fabs (re_tau+0.5);
	return dist;
#endif

}

double domain (double re_q, double im_q)
{
	/* get tau */
	double phase = atan2 (im_q, re_q);
	double r = sqrt (re_q*re_q + im_q*im_q);
	if (1.0<=r) return 0.5;

	r = log (r);
	
	double re_tau = phase / (2.0*M_PI);
	double im_tau = -r/(2.0*M_PI);

	double rc = domain_bounce (re_tau, im_tau, bound, 3);

	return rc;
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
			double re_c = re_position;
			double im_c = im_position;

// #define Q_SERIES_MOBIUS
#ifdef Q_SERIES_MOBIUS
			/* First, make a map from q-series coords to the 
			 * upper half-plane, then apply the mobius x-form, 
			 * and then go back to the q-series coords */
			double qre = re_c;
			double qim = im_c;
			double tau_im = -log (sqrt (qre*qre +qim*qim)) / (2.0*M_PI);
			double tau_re = atan2 (qim, qre) /(2.0*M_PI);

			/* now apply mobius */
			mobius_xform (1, 0, 1, 1, tau_re, tau_im, &tau_re, &tau_im);
			// mobius_xform (0, -1, 1, 0, tau_re, tau_im, &tau_re, &tau_im);

			plane_to_qdisk_coords (tau_re, tau_im, &re_c, &im_c);

#endif /* Q_SERIES_MOBIUS */

			double tau_re, tau_im;
			poincare_disk_to_plane_coords (re_c, im_c, &tau_re, &tau_im);

			mobius_xform (1, 0, 1, 1, tau_re, tau_im, &tau_re, &tau_im);
			// mobius_xform (1, 1, 0, 1, tau_re, tau_im, &tau_re, &tau_im);
			// mobius_xform (0, -1, 1, 0, tau_re, tau_im, &tau_re, &tau_im);
double phi=1.0;
if (tau_re < -0.5) phi = 0.0;
if (tau_re > 0.5) phi = 0.0;
if (tau_re*tau_re+tau_im*tau_im < 1.0) phi=0.0;

			plane_to_q_disk_coords (tau_re, tau_im, &re_c, &im_c);

			// double phi = erdos_series (re_c, im_c);
			// double phi = gee_2 (re_c, im_c);
			// double phi = gee_3 (re_c, im_c);
			// double phi = discriminant (re_c, im_c);
			// double phi = klein_j (re_c, im_c);
			// double phi = domain (re_c, im_c);
         glob [i*sizex +j] = phi;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */
