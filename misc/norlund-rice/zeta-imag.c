/*
 * zeta-imag.c
 *
 * grap[h of zeta on the critical line
 *
 * Linas Vepstas, February 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "harmonic.h"

int
main (int argc, char * argv[])
{
	double ims;
	for (ims=0.0; ims<45; ims +=0.1)
	{
		double rez, imz;
		riemann_zeta (0.5, ims, &rez, &imz);

		printf ("%g\t%g\t%g\n", ims, rez, imz);
	}

	return 0;
}
