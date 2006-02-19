
/*
 * baez.c
 *
 * Baez-Duarte sum
 *
 * Linas Vepstas, February 2006
 */

#include <gsl/gsl_sf_gamma.h>

void integrand (double t, int n, double *reg, double * img)
{
}

main ()
{
	double t=0.0;
	int n=3;

	for (t=-10.0; t<10.0; t+=0.3)
	{
		double reg, img;
		integrand (t, n, &reg, &img);

		printf ("duude its %g  - -%g %g \n", t, reg, img);
	}
}
