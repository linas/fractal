
/*
 * baez.c
 *
 * Baez-Duarte sum
 *
 * Linas Vepstas, February 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_gamma.h>

#include "binomial.h"
#include "harmonic.h"

void log_pole_integrand (double res, double ims, double *reg, double * img)
{
	double regr = 0.0;
	double imgr = 0.0;

	double r = res*res + ims*ims;
	regr -= 0.5*log (r);
	imgr -= atan2 (ims, res);

	/* zeta (2s+2) */
	double rez, imz;
	riemann_zeta (2.0*res+2.0, 2.0*ims, &rez, &imz);

	r = rez*rez + imz*imz;	
	regr -= 0.5 * log (r);
	imgr -= atan2 (imz, rez);

	*reg = regr;
	*img = imgr;
}

void pole_integrand (double res, double ims, double *reg, double * img)
{
	double regr = 1.0;
	double imgr = 0.0;

	double r = res*res + ims*ims;
	double rer = res / r;
	double imr = -ims / r;
	double tmp = regr * rer - imgr * imr;
	imgr = regr * imr + imgr * rer;
	regr = tmp;

	/* zeta (2s+2) */
	double rez, imz;
	riemann_zeta (2.0*res+2.0, 2.0*ims, &rez, &imz);

	r = rez*rez + imz*imz;	

	rer = rez / r;
	imr = -imz / r;
	tmp = regr * rer - imgr * imr;
	imgr = regr * imr + imgr * rer;
	regr = tmp;

	*reg = regr;
	*img = imgr;
}

/* Sanity check the cauchy integral */
double
cauchy_integral (double re_center, double im_center, double radius)
{
	double reacc = 0.0;
	double imacc = 0.0;

	double step = 0.01;
	double theta;
	double npts = 0.0;
	for (theta=-M_PI; theta<=M_PI; theta += step)
	{
		double res = re_center + radius * cos (theta);
		double ims = im_center + radius * sin (theta);

		double reg, img;
		pole_integrand (res, ims, &reg, &img);
		// printf ("%g\t%g\t%g\n", theta, reg, img);

		/* Multiply by the line-element, which is 
		 * tangent to the radial point 
		 */
		double redl = -ims;
		double imdl = res;
		reacc += reg*redl-img*imdl;
		imacc += img*redl+reg*imdl;

		npts += 1.0;
		// printf ("%g\t%g\t%g\n", t, step*reacc, step*imacc);
	}

	reacc /= npts;
	imacc /= npts;

	double tmp = reacc;
	reacc = imacc;
	imacc = -tmp;
	
	printf ("# duude cauchy integral=%g  %g\n", reacc, imacc);
	return reacc;
}

void simple_integrand (double res, double ims, int n, double *reg, double * img)
{
	gsl_sf_result lnr, arg;

	double regr = 0.0;
	double imgr = 0.0;

	/* Gamma (s) */
	gsl_sf_lngamma_complex_e (res, ims, &lnr, &arg);
	regr += lnr.val;
	imgr += arg.val;

	/* Gamma (s+n+1) */
	gsl_sf_lngamma_complex_e (res+n+1.0, ims, &lnr, &arg);
	regr -= lnr.val;
	imgr -= arg.val;

#ifdef VALIDATE
	double ra = exp (regr);
	double x = ra*cos (imgr);
	double y = ra*sin (imgr);
	double f = res*x - ims*y;
	double o = res*y + ims*x;
	printf ("duude s=(%g %g) gam=%g %g\n", res, ims, f, o);
#endif

	/* zeta (2s+2) */
	double rez, imz;
	riemann_zeta (2.0*res+2.0, 2.0*ims, &rez, &imz);
	double r = rez*rez + imz*imz;	

	regr -= 0.5 * log (r);
	imgr -= atan2 (imz, rez);

	*reg = regr;
	*img = imgr;
}

void reflect_integrand (double t, int n, double *reg, double * img)
{
	gsl_sf_result lnr, arg;

	double regr = 0.0;
	double imgr = 0.0;

	double res = 0.5;
	double ims = t;

	/* Gamma (s) */
	gsl_sf_lngamma_complex_e (res, ims, &lnr, &arg);
	regr += lnr.val;
	imgr += arg.val;

	/* Gamma (s+n+1) */
	gsl_sf_lngamma_complex_e (res+n+1.0, ims, &lnr, &arg);
	regr -= lnr.val;
	imgr -= arg.val;

	/* Gamma (2*s+1) */
	gsl_sf_lngamma_complex_e (2.0*res+1.0, 2.0*ims, &lnr, &arg);
	regr -= lnr.val;
	imgr -= arg.val;

	/* (2pi)^(2s-1) */
	double lp = log (2.0*M_PI);
	regr += (2.0*res-1.0)*lp;
	imgr += 2.0*ims*lp;

	/* sin (pi s) */
	double r = exp (-M_PI*ims);
	double rr = 1.0/r;
	double x = (r - rr) * cos (M_PI*res);
	double y = (r + rr) * sin (M_PI*res);

	double theta = atan2 (y, x);
	r = 0.5 * log (x*x+y*y);

	regr -= r;
	imgr -= theta;

	regr += M_LN2;

	/* zeta (2s+1) */
	double rez, imz;
	riemann_zeta (2.0*res+1.0, 2.0*ims, &rez, &imz);
	r = rez*rez + imz*imz;	

	regr -= 0.5 * log (r);
	imgr -= atan2 (imz, rez);

	*reg = regr;
	*img = imgr;
}

double
arc_integral (int n, double re_center, double im_center, double radius)
{
	double reacc = 0.0;
	double imacc = 0.0;

	double step = 0.08;
	double theta;
	double npts = 0.0;
	for (theta=-M_PI; theta<=M_PI; theta+=step)
	{
		double res = re_center + radius * cos (theta);
		double ims = im_center + radius * sin (theta);

		double reg, img;
		// simple_integrand (res, ims, n, &reg, &img);
		pole_integrand (res, ims, &reg, &img);
		// printf ("%g\t%g\t%g\n", t, reg, img);

		double r = exp (reg);
		reacc += r * cos (img);
		imacc += r * sin (img);

		npts += 1.0;
		// printf ("%g\t%g\t%g\n", t, step*reacc, step*imacc);
	}

	reacc /= npts;
	reacc *= factorial (n);
	reacc /= 2.0*M_PI;

	imacc /= npts;
	imacc *= factorial (n);
	imacc /= 2.0*M_PI;

	printf ("# duude arc integral=%g  %g\n", reacc, imacc);
	return reacc;
}

double
integrate (int n)
{
	double t=0.0;

	double reacc = 0.0;
	double imacc = 0.0;

	double step = 0.1;
	double lim = 10.0;
	double npts = 0.0;
	for (t=-lim; t<=lim; t+=step)
	{
		double reg, img;
		simple_integrand (-0.5, t, n, &reg, &img);
		// printf ("%g\t%g\t%g\n", t, reg, img);

		double r = exp (reg);
		reacc += r * cos (img);
		imacc += r * sin (img);

		// printf ("%g\t%g\t%g\n", t, step*reacc, step*imacc);
		npts += 1.0;
	}

	reacc /= npts;
	reacc *= factorial (n);
	reacc /= 2.0*M_PI;

	imacc /= npts;
	imacc *= factorial (n);
	imacc /= 2.0*M_PI;

	printf ("# duude line integral=%g  %g\n", reacc, imacc);
	return reacc;
}

double sum (int n)
{
	int k;
	double acc = 0.0;

	double sgn = 1.0;
	for (k=0; k<=n; k++)
	{
		double b = binomial (n,k);
		double z = zetam1 (2*k+2);
		z = 1.0 / (1.0+z);
		acc += sgn*b*z;
		sgn = -sgn;
	}

	return acc;
}

int
main (int argc, char * argv[])
{
	int n=3;

	if (2>argc)
	{
		fprintf (stderr, "Usage: %s <parm>\n", argv[0]);
		exit (1);
	}
	double rad = atof (argv[1]);
n=0;;

	double in = integrate (n);
	// double ain = arc_integral (n, 0.0, 0.0, rad);
	double ain = cauchy_integral (0.0, 0.0, rad);
	double su = sum (n);

	printf ("# integ=%g arc=%g  sum=%g r = %g\n", in, ain, su, su/in);

#if 0
	for (n=1; n<40; n++)
	{
		su = sum (n);

		double norm = pow (n, 0.75);
		norm *= log(n)*log(n);
		su *= norm;
		printf ("%d\t%g\n", n, su);
	}
#endif

	return 0;
}
