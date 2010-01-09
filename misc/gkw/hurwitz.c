
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
	double s=2.0;

	s = atof (argv[1]);

	int npts = 300;
	double norm = 1.0 / gsl_sf_hzeta(s, 1.0);
	for (i=0; i< npts; i++)
	{
		double x = ((double) i) / ((double) npts);
		double y = norm * gsl_sf_hzeta(s, 1.0+x);

		printf("%d	%g	%g\n", i, x, y);
	}
}
