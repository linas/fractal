
/*
 * dsubn.c
 *
 * proide analytic fit to d_n
 *
 * Linas November 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_sf_gamma.h>
#include "harmonic.h"

void doit(double res, double ims)
{
	gsl_sf_result lnr, arg;
	gsl_sf_lngamma_complex_e (-res, -ims, &lnr, &arg);

	double r = exp(lnr.val);
	printf("Gamma(%g+i%g)= %15.10g * exp(i %15.10g)\n", -res, -ims, r, arg.val);

	double rez, imz, rezz, imzz;
	double h=1.0e-8;
	riemann_zeta (res+h, ims, &rez, &imz);
	riemann_zeta (res, ims, &rezz, &imzz);
	double zpre = (rez-rezz)/h;
	double zpim = (imz-imzz)/h;
	printf ("zeta-prime(%g+i%g)= %15.10g +i %15.10g\n", res, ims, zpre, zpim);

	double zpm = sqrt(zpre*zpre + zpim*zpim);
	double zarg = atan2 (zpim, zpre);

	printf ("zeta-prime(%g+i%g)= %15.10g *exp(i %15.10g)\n", res, ims, zpm, zarg);
	
	double trm = r/zpm;
	double trph = arg.val - zarg;

	printf ("term(%g+i%g)= %15.10g *exp(i %15.10g)\n", res, ims, 2.0*trm, trph);
	
}
	
int
main ()
{

	double res=0.5;
	double ims = 14.13472514173469379;
	
	doit (res, ims);
	
	ims = 21.02203963877;

	doit (res, ims);

	return 0;
}
