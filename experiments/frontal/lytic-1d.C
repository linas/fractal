/*
 *
 * Numerical exploration of the analytic continuation
 * given in frontal.lyx
 *
 * September 2015
 */

#include <stdio.h>
#include <stdlib.h>

#include <complex.h>

#include "lytic.h"

int main(int argc, char *argv[])
{
	double x = 0.49999999999;

	// int k = atoi(argv[1]);
	x = atof(argv[1]);

	printf("#\n# expanding x=%g\n#\n", x);

	for (double z=0.0; z< 1.0; z+= 1.0/1621.0)
	{
		double complex sum = sum_extend(x, z-x);
		printf("%g\t%g\n", z, creal(sum));
	}

	return 0;
}
