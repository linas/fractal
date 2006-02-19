
/*
 * baez.c
 *
 * Baez-Duarte sum
 *
 * Linas Vepstas, February 2006
 */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_sf_gamma.h>

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

	*reg = regr;
	*img = imgr;
}

int
main (int argc, char * argv[])
{
	double t=0.0;
	int n=3;

	for (t=-10.0; t<10.0; t+=0.3)
	{
		double reg, img;
		integrand (t, n, &reg, &img);

		printf ("duude its %g   %g %g \n", t, reg, img);
	}

	return 0;
}
