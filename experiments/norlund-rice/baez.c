
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

/* ======================================================== */
/* Simple test function -- logarithm of simple pole, residue 1, */
void log_pole_integrand (double res, double ims, double *reg, double * img)
{
	double regr = 0.0;
	double imgr = 0.0;

	double r = res*res + ims*ims;
	regr -= 0.5*log (r);
	imgr -= atan2 (ims, res);

	res -= 1.0;

	r = res*res + ims*ims;
	regr -= 0.5*log (r);
	imgr -= atan2 (ims, res);

#if 0
	/* zeta (2s+2) */
	double rez, imz;
	riemann_zeta (2.0*res+2.0, 2.0*ims, &rez, &imz);

	r = rez*rez + imz*imz;	
	regr -= 0.5 * log (r);
	imgr -= atan2 (imz, rez);
#endif

	*reg = regr;
	*img = imgr;
}

/* ======================================================== */
/* An integrand with one single simple pole 
 * having a residue of 1/zeta(2) 
 * Used a sa sanity check for the integral 
 */
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

/* ======================================================== */
/* Sanity check the cauchy integral of a simple pole */
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

	/* Divide by i; the 2pi already came with the line elt */
	double tmp = reacc;
	reacc = imacc;
	imacc = -tmp;
	
	printf ("# duude cauchy integral=%g  %g\n", reacc, imacc);
	return reacc;
}

/* ======================================================== */
/* Integrand corresponding to naive norlund-rice of 
 * the Baez-duarte sum.
 */
void simple_integrand (double res, double ims, int n, double *reg, double * img)
{
	gsl_sf_result lnr, arg;

	double regr = 0.0;
	double imgr = 0.0;

	/* Gamma (s-n) */
	gsl_sf_lngamma_complex_e (res-n, ims, &lnr, &arg);
	regr += lnr.val;
	imgr += arg.val;

	/* Gamma (s+1) */
	gsl_sf_lngamma_complex_e (res+1.0, ims, &lnr, &arg);
	regr -= lnr.val;
	imgr -= arg.val;

// #define VALIDATE
#ifdef VALIDATE_XXXX
	double ra = exp (regr);
	double x = ra*cos (imgr);
	double y = ra*sin (imgr);
	int k;
	for (k=0; k<=n; k++)
	{
		double f = (res+k)*x - ims*y;
		double o = (res+k)*y + ims*x;
		x = f;
		y = o;
	}
	printf ("duude s=(%g %g) pole gam=%g %g\n", res, ims, x, y);
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

/* ======================================================== */
/* Integrand corresponding to norlund-rice of 
 * the Baez-duarte sum. Here, the functional equation for the
 * Riemann zeta has been applied.
 *
 * The real part appears to agree with the non-reflected integrand.
 * The imaginary part seems to be off by pi/2 
 */
void reflect_integrand (double res, double ims, int n, double *reg, double * img)
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

	/* Gamma (2*s-1) */
	gsl_sf_lngamma_complex_e (2.0*res-1.0, 2.0*ims, &lnr, &arg);
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

	// regr += M_LN2;

	/* zeta (2s-1) */
	double rez, imz;
	riemann_zeta (2.0*res-1.0, 2.0*ims, &rez, &imz);
	r = rez*rez + imz*imz;	

	regr -= 0.5 * log (r);
	imgr -= atan2 (imz, rez);

	*reg = regr;
	*img = imgr;
}

/* ======================================================== */
/* Integrate along a half-circle arc */
double
arc_integral (int n, double re_center, double im_center, double radius)
{
	double reacc = 0.0;
	double imacc = 0.0;

	double step = 0.003;
	double theta;
	double npts = 0.0;
	for (theta=-0.5*M_PI; theta<=0.5*M_PI; theta+=step)
	{
		double res = re_center + radius * cos (theta);
		double ims = im_center + radius * sin (theta);

		double reg, img;
		simple_integrand (res, ims, n, &reg, &img);
		// log_pole_integrand (res, ims, &reg, &img);
		// printf ("%g\t%g\t%g\n", t, reg, img);

		double r = exp (reg);
		double rei = r * cos (img);
		double imi = r * sin (img);

		/* Multiply by the line-element, which is 
		 * tangent to the radial point 
		 */
		double redl = -ims;
		double imdl = res;
		reacc += rei*redl-imi*imdl;
		imacc += imi*redl+rei*imdl;

		npts += 1.0;
		// printf ("%g\t%g\t%g\n", t, step*reacc, step*imacc);
	}

	reacc /= npts;
	reacc *= factorial (n);

	imacc /= npts;
	imacc *= factorial (n);

	/* Divide by 2i; the pi already came with the line elt */
	double tmp = 0.5*reacc;
	reacc = 0.5*imacc;
	imacc = -tmp;

	printf ("# duude arc integral=%g  %g\n", reacc, imacc);
	return reacc;
}

/* ======================================================== */
/* integrate along imaginary direction */

double
integrate (int n, double re_offset, double lim)
{
	double t=0.0;

	double reacc = 0.0;
	double imacc = 0.0;

	double step = 0.002*lim;
	double npts = 0.0;
	double imlast = 0.0;
	double lgf = log (factorial (n));
	lgf -= log (2.0*M_PI);
	for (t=-lim; t<=lim; t+=step)
	{
		double reg, img;
		simple_integrand (re_offset, t, n, &reg, &img);
		// simple_integrand (-0.25, t, n, &reg, &img);

		// double nreg, nimg;
		// reflect_integrand (re_offset, t, n, &reg, &img);
		// reflect_integrand (0.5, t, n, &nreg, &nimg);
		// log_pole_integrand (re_offset, t, &reg, &img);

		while (img >imlast+5.0) img -= 2.0*M_PI;
		while (img <imlast-5.0) img += 2.0*M_PI;
		imlast = img;

		reg += lgf;

		printf ("%g\t%g\t%g\n", t, reg, img);
		// printf ("%g\t%g\t%g   %g\n", t, img, nimg, img+nimg);

		double r = exp (reg);
		double ret = r * cos (img);
		double imt = r * sin (img);
		// printf ("%g\t%g\t%g\n", t, ret, imt);

		reacc += ret;
		imacc += imt;

		// printf ("%g\t%g\t%g\n", t, step*reacc, step*imacc);
		npts += 1.0;
	}

	reacc *= 2.0*lim/npts;
	// reacc *= factorial (n);
	// reacc /= 2.0*M_PI;

	imacc *= 2.0*lim/npts;
	// imacc *= factorial (n);
	// imacc /= 2.0*M_PI;

	printf ("# duude line integral re=%g  im=%g\n", reacc, imacc);
	return reacc;
}

/* ======================================================== */
/* just print data */

void
show_integrand (int n, double im_offset, double lim)
{
	double t=0.0;

	double step = 0.001;
	double imlast = 0.0;
	double lgf = log (factorial (n));
	lgf -= log (2.0*M_PI);

	double prev = 0.0;
	for (t=-0.1; t>-1.9; t-=step)
	{
		double reg, img;
		simple_integrand (t, im_offset, n, &reg, &img);

		while (img >imlast+5.0) img -= 2.0*M_PI;
		while (img <imlast-5.0) img += 2.0*M_PI;
		imlast = img;

		reg += lgf;

		// printf ("%g\t%g\t%g\n", t, reg, img);
		// printf ("%g\t%g\t%g   %g\n", t, img, nimg, img+nimg);

		double r = exp (reg);
		double ret = r * cos (img);
		double imt = r * sin (img);

		printf ("%g\t%g\t%g\n", t, ret, (prev-ret)/step);
		prev = ret;
	}
}

/* ======================================================== */
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

/* ======================================================== */
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
	n=rad;

	double offset = rad;
	n = 26;
	rad = 14.0;

	// double in = integrate (n, offset, rad);
	// double ain = arc_integral (n, -0.5, 0.0, rad);
	// double ain = cauchy_integral (0.0, 0.0, rad);
	show_integrand (n, 0.0, rad);
	double su = sum (n);

	// printf ("# integ=%g arc=%g  sum=%g r = %g\n", in, ain, su, );
	// printf ("# integ=%g sum=%g r = %g\n", in, su, in/su);

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
