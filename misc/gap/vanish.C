
/*
 * vanish.C
 *
 * functions that vanish on binary numbers
 */

#include <math.h>
#include <stdio.h>

double sin_adic (int n, double x)
{
	int i;
	double acc = 1.0;
	double tp = M_PI;
	for (i=0; i<n; i++)
	{
		acc *= sin (tp*x);
		tp *= 2.0;
	}
	return acc;
}

main ()
{
	int i;

	int nmax=523;
	for (i=0; i<nmax; i++)
	{
		double x = i/((double) nmax);

		double y = sin_adic (20, x);
		printf ("%d	%8.6g	%8.6g\n", i, x, y);
	}
}
