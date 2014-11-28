//
// bessy.cc
//
// Solve schroedinger eqn in a 1D lattice.  Graph the bessel-function-like
// solutions.
//
// Linas Vepstas 23 Nov 2014
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "moebius.h"

#include "schroed.h"

// Well, if its bessel-like, then lets go with that idea...
void bessy(size_t len, double omega, double maxen)
{
	long double *p = (long double*) malloc (len * sizeof(long double));
	for (size_t i=0; i<len; i++) p[i] = sigmaf(i, omega);
	set_pot (len, p);

	// printf("#\n#Bessel-like totient solution\n#\n");
	printf("#\n# Divisor solution\n#\n");
	printf("# Energy  discontinuity b0	b1	b2	b3	b4	b5\n");
	double en=0.0;
	double delta = maxen / 1000.0;
	for (en=0.0; en<maxen; en+=delta)
	{
		long double discon = solve(len, 1.0, en);

		printf("%g	%Lg	%Lg	%Lg	%Lg	%Lg	%Lg	%Lg\n",
			en, discon, wavefn[0], wavefn[1], wavefn[2], wavefn[3],
			wavefn[4], wavefn[5]);
	}
}

main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s <num-grid-pts> <max-energy> <omega>\n", argv[0]);
		exit (1);
	}
	int len = atoi(argv[1]);

	double maxen = atof(argv[2]);
	double omega = atof(argv[3]);

	init(len, omega);

   // solve(len, maxen);
	bessy(len, omega, maxen);
}
