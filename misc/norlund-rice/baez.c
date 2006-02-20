
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
#include <gsl/gsl_sf_zeta.h>

#include "binomial.h"

/* Compute the Riemann zeta for general complex argument.
 * Uses the Hasse expansion. More or less accurate;
 * has some accuracy trouble near the pole at s=1;
 */
void riemann_zeta (double res, double ims, double *rez, double *imz)
{
	double zre = 0.0;
	double zim = 0.0;
	double err;

	int n;

	double twn = 0.5;

	/* don't change the 110 -- it seems liek the right 
	 * thing for stuff going in the imaginary direction
	 */
	for (n=0; n<110; n++)
	{
		int k;
		double reb = 0.0;
		double imb = 0.0;
		double sgn = 1.0;
		for (k=0; k<=n; k++)
		{
			double r = binomial (n,k);
			// printf ("duude %d %d bin=%g\n", n, k, r);
			double lnk = log (k+1.0);
			r *= sgn * exp (-res*lnk);
			reb += r * cos (ims*lnk);
			imb -= r * sin (ims*lnk);
			sgn = -sgn;
		}

		double ret = twn * reb;
		double imt = twn * imb;
		zre += ret;
		zim += imt;

		err = ret*ret+imt*imt;
		// printf ("duude n=%d z=(%g %g) err=%g\n", n, zre, zim, err);

		/* Along imaginary axies, the error never seems 
		 * get less than 1.0e-34, no matter what (I don't get why)
		 * so cut off here.
		 */
		if (err < 1.0e-33) break;
		twn *= 0.5;
	}

	// printf ("duude leave with n=%d err=%g \n", n, sqrt (err));

	double r = pow (2.0, 1-res);
	double ret = 1.0 - r* cos (ims*M_LN2);
	double imt = r* sin (ims*M_LN2);
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
	double s,t=0.0;
	int n;

	int error_occured = 0;

	for (n=2; n<=40; n++)
	{
		double reg, img;
		riemann_zeta (n, 0.0, &reg, &img);
		double gslz = gsl_sf_zeta_int (n);

		double err = reg-gslz;
		if ((fabs(err) > 1.0e-15) || (fabs (img) > 1.0e-15))
		{
			printf ("ERROR for n=%d   error=%g %g \n", n, err,  img);
			error_occured ++;
		}
	}

	for (t=2.0; t<=66.0; t+=0.0314683)
	{
		double reg, img;
		riemann_zeta (t, 0.0, &reg, &img);
		double gslz = gsl_sf_zeta (t);

		double err = reg-gslz;
		if ((fabs(err) > 1.0e-15) || (fabs (img) > 1.0e-15))
		{
			printf ("ERROR for s=%g   error=%g %g \n", t, err,  img);
			error_occured ++;
		}
	}

	for (t=1.01; t<=2.1; t+=0.00314683)
	{
		double reg, img;
		riemann_zeta (t, 0.0, &reg, &img);
		double gslz = gsl_sf_zeta (t);

		double err = reg-gslz;
		if ((fabs(err) > 1.0e-13) || (fabs (img) > 1.0e-15))
		{
			printf ("ERROR for s=%g   error=%g %g \n", t, err,  img);
			error_occured ++;
		}
	}

	for (s=-40.0; s<=40.0; s += 0.4356346)
	{
		for (t=0.0; t<=48.0; t+=0.6314683)
		{
			double reg, img;
			riemann_zeta (0.5, t, &reg, &img);
			double nreg, nimg;
			riemann_zeta (0.5, -t, &nreg, &nimg);
	
			double rerr = reg-nreg;
			double ierr = img+nimg;
			if ((fabs(rerr) > 1.0e-13) || (fabs (ierr) > 1.0e-15))
			{
				printf ("ERROR for s=%g   error=%g %g \n", t, rerr,  ierr);
				error_occured ++;
			}
		}
	}

	return error_occured;
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

		printf ("%g\t%g\t%g\n", t, reg, img);
	}

	return 0;
}
#endif
