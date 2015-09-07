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

	for (double u=-1.0; u< 1.0; u+= 1.0/621.0)
	{
		printf("%g", u);
		for (int k=1; k<6; k++)
		{
			double complex cu_k = count_extend(k, x, u);
			printf("\t%g", creal(cu_k));
		}
		printf("\n");
	}

	return 0;
}
