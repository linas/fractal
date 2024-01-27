/*
 * phasor.C
 *
 * Version of unstack.c with complex lambda inserted into it.
 * Explore behavior on circle.
 *
 * January 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// ==============================================================

// Arbitrary function
double nu(double x)
{
	if (x < 0.0) fprintf(stderr, "Error nu fail neg %g\n", x);
	if (1.0 < x) fprintf(stderr, "Error nu fail pos %g\n", x);

	// return 1.0;
	// return x-0.5;
	return x - 0.5 + 0.08684;  // appropriate for beta=1.6

	// Bernoulli poly B_2
	// The result is senstive to this being B_2.
	// Being able to integrate to exactly zero is important.
	// return x*x - x  + 1.0 / 6.0;
	// return x*x - x  + 0.16666;

	// Bernoulli poly B_3
	// return x*x*x - 1.5*x*x  + 0.5*x;

	// Bernoulli poly B_4
	// return x*x*x*x - 2.0*x*x*x  + x*x - 1.0/30.0;
}

#include "uncomplex.C"
#include "unlambda.c"
#include "unref.c"

// ==============================================================

int main(int argc, char* argv[])
{
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s beta delta-phi n\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	double dphi = atof(argv[2]);
	int n = atoi(argv[3]);

	// double scale = lambda * beta;
	// scale = lambda;
	double scale = 1.0;
	double scan = pow(scale, n);
	int imax = 814;
	for (int i=0; i< imax; i++)
	{
		double x = (((double) i) + 0.5) / ((double) imax);

		double y = gp_invar(beta, x);
		printf("%d	%f	%f", i, x, y);

		double plm = scan;
#define NIT 6
		for (int j=0; j<NIT; j++)
		{
			COMPLEX blam = cexp(2.0*M_PI*I*j*dphi);
			COMPLEX y = nuz_n(beta, blam, x, n);
			// y /= blam;
			// double why = nul_n(beta, 1.0, x, n);
			// printf ("	%f", why);
			y *= plm;
			plm *= scale;
			printf("	%f", std::real(y));
			printf("	%f", imag(y));
		}
		printf("\n");
		fflush(stdout);
	}
}
