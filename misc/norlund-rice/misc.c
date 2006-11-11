
#include <gsl/gsl_sf_gamma.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
	gsl_sf_result lnr, arg;
	
	gsl_sf_lngamma_complex_e (-0.5, 14.13472514173469379, &lnr, &arg);

	double r = exp(lnr.val);

	printf("duude its %15.10g * exp(i %15.10g)\n", r, arg.val);

	gsl_sf_lngamma_complex_e (-0.5, 21.02203963877, &lnr, &arg);

	r = exp(lnr.val);

	printf("duude its %15.10g * exp(i %15.10g)\n", r, arg.val);

	return 0;
}
