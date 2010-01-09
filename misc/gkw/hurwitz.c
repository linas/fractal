
/*
 * graph of hurwitz zeta
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_sf_zeta.h>


main ()
{
	double s=2.o;

	int npts = 300;
	for (int i=0; i< npts; i++)
	{
		double x = ((double) i) / ((double) npts);
		double y = gsl_sf_hzeta(2.0, x)

		print ("%d	%g	%g\n", i, x, y);
	}
}
