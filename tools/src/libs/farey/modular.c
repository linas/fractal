/*
 * modular.C
 *
 * FUNCTION:
 * assorted modular functions.
 * -- divisor arithmetic function
 * -- Weierstrass elliptic function g_2 and g_3 invariants
 * -- modular discriminant.
 * -- Klein j-invariant
 *
 * HISTORY:
 * more stuff -- October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "modular.h"
#include "moebius.h"

int modular_max_terms;

/** Compute the divisor arithmetic function
 *  Returns the number of divisors of n.
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

/** Sigma arithmetic series, equals divisor airth series for a=0 
 *  Computes the divisors of n, raises each to the a'th power, and
 *  returns thier sum.
 */
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

/* An erdos-borwein-like series sum_n sigma_k(n) x^n/(1-x^n) 
 * Actually computed via a q-series for the divisor function.
 * This is a building block for modular forms.
 */
void erdos_series_c (long double re_q, 
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
	for (i=0; i<modular_max_terms; i++)
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
	if (modular_max_terms-1 < i)
	{
		// printf ("not converged re=%g im=%g modulus=%g\n", re_q, im_q, qpmod);
	}

	*prep = rep;
	*pimp = imp;
}

/* Weierstrass elliptic invarient g_2, where q is the nome */
void gee_2_c (long double re_q, long double im_q, long double *pre, long double *pim)
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

void gee_3_c (long double re_q, long double im_q, long double *pre, long double *pim)
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

/** The modular discriminant, constructed as g_2^3 - 27g_3^2 . 
 *  Not very accurate, due to the need to take teh difference of two
 *  paritally converged sums.
 */
void disc_c (long double re_q, long double im_q, long double *pre, long double *pim)
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

	*pre = dre;
	*pim = dim;
}

/** The modular discriminat, computed as dedekind eta to 24 
 */
void discriminant_c (double re_q, double im_q, double *pre,double *pim)
{
	double rep, imp;
	euler_prod_c (re_q, im_q, &rep, &imp);

	/* Take euler product to 24'th power */
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

	*pre = red;
	*pim = imd;
}

/* Return euler-product form of the q-series (dedekind eta) */
void euler_prod_c (double re_q, double im_q, double *prep, double *pimp)
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

/* The dedekind eta multiplies by an additonal factor of q^1/24 */
void dedekind_eta_c (double re_q, double im_q, double *pre, double *pim)
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

/* The Klein j invariant, */
void klein_j_invariant_c (double re_q, double im_q, double *pre, double *pim)
{
	double tmp, dre, dim;

	/* compute 1/delta */
	discriminant_c (re_q, im_q, &dre, &dim);
	double dmod = dre*dre+dim*dim;
	dre /= dmod;
	dim = -dim/dmod;

	/* compute g_2 cubed */
	long double g2re, g2im;
	gee_2_c (re_q, im_q, &g2re, &g2im);
	
	long double g2cure, g2cuim;
	g2cure = g2re*g2re - g2im*g2im;
	g2cuim = 2.0L * g2re*g2im;
	tmp = g2cure*g2re - g2cuim*g2im;
	g2cuim = g2cure*g2im + g2cuim* g2re;
	g2cure = tmp;

	/* compute g_2 cubed /delta */
	double jre = dre*g2cure - dim*g2cuim;
	double jim = dre*g2cuim + dim*g2cure;

	*pre = jre;
	*pim = jim;
}

/* --------------------------- END OF FILE ------------------------- */
