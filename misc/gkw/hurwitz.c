
/*
 * graph of hurwitz zeta
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_sf_zeta.h>


main (int argc, char * argv[])
{
	int i;
	double s = 2.0;

	s = atof (argv[1]);

	int npts = 300;
	double n2 = gsl_sf_hzeta(2.0, 1.0);
	double n3 = gsl_sf_hzeta(3.0, 1.0);
	double n4 = gsl_sf_hzeta(4.0, 1.0);

	// double norm = 1.0 / (n2 - s * n3);
	double norm = 1.0 / (n2 - s * n4);
	for (i=0; i< npts; i++)
	{
		double x = ((double) i) / ((double) npts);
		double y2 = gsl_sf_hzeta(2.0, 1.0+x);
		// double y3 = gsl_sf_hzeta(3.0, 1.0+x);
		double y4 = gsl_sf_hzeta(4.0, 1.0+x);
		// double y = norm * (y2 - s * y3);
		double y = norm * (y2 - s * y4);

		printf("%d	%g	%g\n", i, x, y);
	}
}
