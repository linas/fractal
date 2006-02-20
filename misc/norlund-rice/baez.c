
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

#include "harmonic.h"

void integrand (double t, int n, double *reg, double * img)
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
	double x = r * cos (M_PI*res);
	double y = r * sin (M_PI*res);

	r = 1.0/r;
	x -= r * cos (M_PI*res);
	y += r * sin (M_PI*res);

	double theta = atan2 (y, x);
	r = 0.5 * log (x*x+y*y);

	regr -= r;
	imgr -= theta;

	/* zeta (2s+1) */
	double rez, imz;
	riemann_zeta (2.0*res+1.0, 2.0*ims, &rez, &imz);
	r = rez*rez + imz*imz;	

	regr -= 0.5 * log (r);
	imgr -= atan2 (imz, rez);

	*reg = regr;
	*img = imgr;
}

int
main (int argc, char * argv[])
{
	double t=0.0;
	int n=3;

	for (t=-10.0; t<=10.0; t+=0.2)
	{
		double reg, img;
		integrand (t, n, &reg, &img);

		printf ("%g\t%g\t%g\n", t, reg, img);
	}

	return 0;
}
