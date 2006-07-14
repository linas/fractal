/*
 * mobius.C
 *
 * FUNCTION:
 * display mobius mu, Euler's totient, and other number-theoretic
 * functions on the complex disk (poincare disk).
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
#include "totient.h"

int thue_morse(int n)
{
	if (0 == n) return 0;
	if (1 == n) return 1;
	if (0 == n%2) return thue_morse (n/2);
	return (1-thue_morse ((n-1)/2));
}


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
			array[i] = (rand()>>6) & 0x1;
			// array[i] = ((double) rand()) / ((double)RAND_MAX);
		}
	}

	if (NVAL<= n) return 0;
	return array[n];
}

static int max_terms;

/* Perform ordinary series sum over one of the funcs */
void plain_series_c (double re_q, double im_q, double *prep, double *pimp)
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
#if 0
		t *= (i+1);
		t *= (i+1);
		t *= (i+1);
#endif
		t *= (i+1);
		t *= (i+1);
		t *= (i+1);
		t *= (i+1);
		t *= (i+1);
		t *= (i+1);
		t *= (i+1);
		t *= (i+1);
		t *= (i+1);

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

/* derivatives */
struct deriv_s 
{
	double reh;
	double imh;
	double rehp;
	double imhp;
	double rehpp;
	double imhpp;
};

static void derivatives (double re_q, double im_q, struct deriv_s *dd)
{
	double tmp;
	int i;

	double reh = 0.0;
	double imh = 0.0;
	double rehp = 0.0;
	double imhp = 0.0;
	double rehpp = 0.0;
	double imhpp = 0.0;

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
		double t = totient_phi (i+1);
		// double t = randoid (i+1);
#if 0
		t *= (i+1);
		t *= (i+1);
		t *= (i+1);
#endif

		double eye = i;

		reh += t*qpr;
		imh += t*qpi;

		rehp += eye*qprm1*t;
		imhp += eye*qpim1*t;
		rehpp += eye*(eye-1.0)*qprm2*t;
		imhpp += eye*(eye-1.0)*qpim2*t;

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

	dd->rehpp = rehpp;
	dd->imhpp = imhpp;
	dd->rehp = rehp;
	dd->imhp = imhp;
	dd->reh = reh;
	dd->imh = imh;
}


/* Computes curvature of geodesics, i.e. curvature of field lines (rays)
 * and of equipotentials.  Curvature of field line returned in second,
 * cuvature of equipotential in first. */
void line_curvature_c (double re_q, double im_q, double *pequi, double *pfie)
{
	*pequi = 0.0;
	*pfie = 0.0;

	struct deriv_s dd;
	derivatives (re_q, im_q, &dd);

	double norm = pow (dd.rehp*dd.rehp+dd.imhp*dd.imhp,  1.5);
	double  equipot = - dd.rehpp*(dd.rehp*dd.rehp - dd.imhp*dd.imhp) - 2.0*dd.rehp*dd.imhp*dd.imhpp;
	equipot /= norm;

	double ray = dd.imhpp*(dd.rehp*dd.rehp - dd.imhp*dd.imhp) - 2.0*dd.rehp*dd.imhp*dd.rehpp;
	ray /= norm;

	*pequi = equipot;
	*pfie = ray;
}


/* Computes scalar curvature of surface (contraction of ricci curvature) 
 */
void scalar_c (double re_q, double im_q, double *pcurv, double *pxxx)
{
	struct deriv_s dd;
	derivatives (re_q, im_q, &dd);

	double deno = 1.0+dd.rehp*dd.rehp+dd.imhp*dd.imhp;
	deno *= deno;

	double  numer = - dd.rehpp*dd.rehpp - dd.imhpp*dd.imhpp;
	double curvature = -2.0*numer / deno;

	*pcurv = curvature;
	*pxxx = 0.0;
}

/* Someday,this is going to Compute components of energy-momentum tensor 
   (actually, just the mass)  right now its a test function.
 */
void energy_c (double re_q, double im_q, double *energy, double *moment)
{
	struct deriv_s dd;
	derivatives (re_q, im_q, &dd);

	double gxx = 1.0+dd.rehp*dd.rehp;
	double gxy = -dd.rehp*dd.imhp;
	double gyy = 1.0+dd.imhp*dd.imhp;

	double deno = 1.0+dd.rehp*dd.rehp+dd.imhp*dd.imhp;

	double flub = dd.imhp*dd.imhp * dd.rehp*dd.rehp;
	flub /= deno;


	*energy = flub;
}


static double 
mobius_series (double re_q, double im_q, int itermax, double param)
{
	double rep, imp;

	max_terms = itermax;

	// plain_series_c (re_q, im_q, &rep, &imp);
	// line_curvature_c (re_q, im_q, &rep, &imp);
	scalar_c (re_q, im_q, &rep, &imp);
	// return sqrt (rep*rep+imp*imp);
	return rep;
	// return rep*imp;
	// return imp;
	// return (atan2 (imp,rep)+M_PI)/(2.0*M_PI);
}

DECL_MAKE_HEIGHT(mobius_series)

/* --------------------------- END OF LIFE ------------------------- */
