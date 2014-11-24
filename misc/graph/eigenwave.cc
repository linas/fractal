//
// eigenwave.cc
//
// Solve schroedinger eqn in a 1D lattice.
//
// Linas Vepstas 23 Nov 2014
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "schroed.h"

// Find zero
long double find_eig(size_t len, long double omega, long double lo, long double hi)
{
	long double dilo = solve(len, omega, lo);
	dilo = wavefn[0];
	long double dihi = solve(len, omega, hi);
	dihi = wavefn[0];

	if (fabsl(dilo) < 1.0e-17) return lo;
	if (fabsl(dihi) < 1.0e-17) return hi;

	// compute the midpoint
	long double slope = (dihi - dilo) / (hi - lo);

	long double xmid = lo - (dilo / slope);
	long double ymid = solve(len, omega, xmid);
	ymid = wavefn[0];

	printf("duuude lo %Lg %Lg\n", lo, dilo);
	printf("duuude mid %Lg %Lg\n", xmid, ymid);
	printf("duuude hi %Lg %Lg\n\n", hi, dihi);
	if (ymid < 0.0L and dilo < 0.0L) return find_eig(len, omega, xmid, hi);
	if (ymid > 0.0L and dilo > 0.0L) return find_eig(len, omega, xmid, hi);

	return find_eig(len, omega, lo, xmid);
}

// Graph of wave function at energy
void graf(size_t len, double omega, long double en)
{
	en = find_eig(len, omega, en-0.1, en+0.1);


	long double discon = solve(len, omega, en);
	printf("#\n# Totient wavefunction at energy=%22.18Lf discon=%Lg\n#\n", en, discon);

	for (size_t i=0; i<len; i++)
	{
		printf("%zd	%Lg\n", i, wavefn[i]);
	}
}

main(int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <num-grid-pts> <energy>\n", argv[0]);
		exit (1);
	}
	int len = atoi(argv[1]);

	double en = atof(argv[2]);
	double omega = 1.0;

	init(len, omega);

	graf(len, omega, en);
}
