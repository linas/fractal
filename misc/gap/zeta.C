
/* 
 * zeta.C
 * 
 * explore Hurwitz zeta eigenfunctions of Bernoulli map
 *
 * Linas November 2004
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double beta (double x, double s)
{
	int n;
	double acc = 0.0;
	for (n=1; n<57550; n++)
	{
		double term = pow(2.0*M_PI*((double) n), -s);
		acc += term * cos (2.0*M_PI*((double) n)*x);
		// acc += term * sin (2.0*M_PI*((double) n)*x);
		if (1.0e-16>term) break;
	}
	acc *= 2.0 * tgamma (s+1.0);
	return acc;
}

main (int argc, char * argv[])
{
	int i;

	double s=3.345;

	if (2>argc)
	{
		printf ("Usage: %s  <s-value>\n", argv[0]);
		exit (1);
	}
	s = atof (argv[1]);

	double lambda = pow (0.5, s);

	printf ("#\n# ess=%g  eigenvalue lambda=%g\n#\n", s, lambda);
	
	int imax = 23;
	for (i=0; i<imax; i++) 
	{
		double x = i/((double) imax);
		double y = beta (x,s);

// #define CHECK_THAT_ZETA_IS_EIGENSTATE 1
#if CHECK_THAT_ZETA_IS_EIGENSTATE 
		// verify its an eigenstate
		double z = 0.5*(beta (0.5*x, s) + beta (0.5+0.5*x, s));
		z /= lambda;
		z -= y;
#endif
		double z=1.0;

		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x, y, z);
	}
}
