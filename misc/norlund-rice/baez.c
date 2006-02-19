
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

/* uses the Hasse expansion */
void riemann_zeta (double res, double ims, double *rez, double *imz)
{
	double zre = 0.0;
	double zim = 0.0;

	int n;

	double twn = 0.5;
	for (n=0; n<50; n++)
	{
		int k;
		double reb = 0.0;
		double imb = 0.0;
		double sgn = 1.0;
		for (k=0; k<=n; k++)
		{
			double lnk = log (k+1.0);
			double r = exp (-res*lnk);
			r *= binomial (n,k);
			r *= sgn;
			reb += r * cos (ims*lnk);
			imb += r * sin (ims*lnk);
			sgn = -sgn;
		}

		double ret = twn * reb;
		double imt = twn * imb;
		zre += ret;
		zim += imt;

		if (ret*ret+imt*imt < 1.0e-30) break;
		twn *= 0.5;
	}

printf ("duude leave with n=%d\n", n);

	double r = pow (2.0, 1-res);
	double ret = 1.0 - r* cos (ims*M_LN2);
	double imt = r* sinc (ims*M_LN2);
	r = ret*ret + imt*imt;
	r = 1.0/r;

	ret *= r;
	imt = -imt*r; 

	double tmp = ret * zre - imt *zim;
	zim = ret*zim + imt *zre;
	zre = tmp;
	
	*rez = zre;
	*imz = zim;
}
int
main (int argc, char * argv[])
{
	double t=0.0;
	int n=3;

	for (t=-10.0; t<=10.0; t+=0.2)
	{
		double reg, img;
		// integrand (t, n, &reg, &img);
		reimann_zeta (t, n, &reg, &img);

		printf ("duude its %g   %g %g \n", t, reg, img);
	}

	return 0;
}

#if 0
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
		// integrand (t, n, &reg, &img);
		reimann_zeta (t, n, &reg, &img);

		printf ("duude its %g   %g %g \n", t, reg, img);
	}

	return 0;
}
#endif
