
/*
 * alpha.c
 *
 * check analytic properties of 
 * alpha(z) = sum_n a_n z^n
 *
 * The result of this is to confirm that alpha seems to be entire,
 * and free of any singularities whatsoever.
 *
 * Linas Vepstas
 * December 2004
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>

#include "ache.h"
#include "zetafn.h"

long double complex alpha_series (long double complex z)
{
	long double complex zn, acc;
	acc = 0.0;
	zn = 1.0;

	int n;
	for (n=0; n<40; n++)
	{
		long double an = a_sub_n (n);
		acc += an * zn;
		zn *= z;
	}
	
	return acc;
}

/* Log[Gamma(z)] for z complex, z not a negative integer
 * Uses complex Lanczos method. Note that the phase part (arg)
 * is not well-determined when |z| is very large, due
 * to inevitable roundoff in restricting to (-Pi,Pi].
 * This will raise the GSL_ELOSS exception when it occurs.
 * The absolute value part (lnr), however, never suffers.
 *
 * Calculates:
 *   lnr = log|Gamma(z)|
 *   arg = arg(Gamma(z))  in (-Pi, Pi]
 *
 * exceptions: GSL_EDOM, GSL_ELOSS
int gsl_sf_lngamma_complex_e(double zr, double zi, gsl_sf_result * lnr, gsl_sf_result * arg);
 */


long double complex alpha_func (long double complex z)
{
	gsl_sf_result lnr, arg;
	long double complex alph;

	long double complex vee;

	vee = 1.0 / (1.0-z);
	alph = clogl (1.0-z) / z;
	alph  *= vee - 0.5;
// printf ("duuude                 %Lg  %Lg\n", creall(alph), cimagl(alph));
	alph += vee;

	double vr = creall (vee);
	double vi = cimagl (vee);
	int rc = gsl_sf_lngamma_complex_e (vr, vi, &lnr, &arg);

	long double complex gam;
	// gam = lnr.val + I * (arg.val-2.0*M_PI);
	gam = lnr.val + I * arg.val;

	gam /= z;

// printf ("duuude   alp= %Lg  %Lg    gam=%Lg   %Lg\n", creall(alph), cimagl(alph), creall(gam), cimagl(gam));
	alph += gam;
	
	return alph;
}

int
main (int argc, char ** argv)
{
	long double complex as, af, z;

	// double one = 0.5*log (2.0*M_PI) - 1.0;
	
	z = 0.95;
	
	int imax = 23;
	int i;
	for (i=1; i<imax; i++)
	{
		z = ((double) i)/ ((double) imax);

		z = cexpl (I*(z-0.5)*0.5*M_PI);

		z += 1.0/51.0;
		z += 3.0;

		as = alpha_series (z);
		af = alpha_func (z);

		printf ("z=%Lg + i%Lg   \taseries=%Lg + i %Lg   \tafunc=%Lg + i %Lg   \tdiff=%Lg\n", 
						creall(z), cimagl (z),
						 creall(as), cimagl(as), creall(af), cimagl(af), cabsl(af-as));
	}
	return 0;
}
