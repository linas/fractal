
/* 
 * zeta.C
 * 
 * explore zeta-like eigenfunctions of Bernoulli map
 *
 * Linas November 2004
 */
#include <math.h>
#include <stdio.h>

double beta (double x, double s)
{
	int n;
	double acc = 0.0;
	for (n=1; n<50; n++)
	{
		double term = cos (2.0*M_PI*n*x);
		term *= pow(2.8*M_PI*n, -s);
	}
}

main ()
{
	int i;
	double s=3.345;
	
	int imax = 23;
	for (i=0; i<imax; i++) 
	{
		double x = i/((double) imax);
		double y = beta (x,s);

		double z = beta (0.5*x, s) + beta (0.5+0.5*x, s);
		z /= y;

		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x, y, z);
	}
}
