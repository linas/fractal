
/*
 * geodesic.c
 *
 * Graphs of the lengths and energies of geodesics on the riemann 
 * surface/fundamental domain. For now,the symmetric orbits only.
 */

#include <math.h>
#include "geo-lib.h"


main () 
{
	int i;

	int n = 1400;
	double rho = 1.0;
	double delta = 0.1;
	rho += delta;

	printf ("#\n# geodesic lengths\n#\n");
	printf ("# delta = %g\n#\n", delta);

	for (i=0; i<n; i++)
	{
		double len = geo_length (rho);
		double eng = geo_energy (rho);
		double ex = geo_x (rho);

		printf ("%d	%g	%g	%g\n", i, rho, len, eng);
		rho += delta;
	}
}
