
/* 
 * period.c
 *
 * understand periodic orbit contstraints
 */

#include <math.h>
#include "geo-lib.h"

double invert_radius (double rho)
{
	double en = get_n_of_rho (rho);

	en = (en+1.0)*(en+1.0);
	double ri = (en + 1.0) / (rho*rho - en);
	ri = sqrt (ri);

	return ri;
}

main () 
{
	int i;

	int n = 400;
	double rho = 1.0;
	double delta = 0.01;
	rho += delta;

	for (i=0; i<n; i++)
	{
		double ri = invert_radius (rho);
		printf ("%d	%g	%g\n", i, rho, ri);
		rho += delta;
	}
}
