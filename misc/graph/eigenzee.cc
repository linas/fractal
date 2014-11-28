//
// eigenzee.cc
//
// Graf of eigenvalues of disvisor function, as function of power.
//
// Linas Vepstas 27 Nov 2014
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "moebius.h"

#include "schroed.h"

// Find zero. It should be brakceted by hi, lo
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

#if 1
	printf("# duuude lo %Lg %Lg\n", lo, dilo);
	printf("# duuude mid %Lg %Lg\n", xmid, ymid);
	printf("# duuude hi %Lg %Lg\n#\n", hi, dihi);
#endif
	if (ymid < 0.0L and dilo < 0.0L) return find_eig(len, omega, xmid, hi);
	if (ymid > 0.0L and dilo > 0.0L) return find_eig(len, omega, xmid, hi);

	return find_eig(len, omega, lo, xmid);
}

#if 0
// Graph of wave function at energy
void graf(size_t len, double omega, long double en)
{
	en = find_eig(len, omega, en-0.01, en+0.01);


	long double discon = solve(len, omega, en);
	printf("#\n# Totient wavefunction at energy=%22.18Lf discon=%Lg\n#\n", en, discon);

	for (size_t i=0; i<len; i++)
	{
		printf("%zd	%22.18Lg	%22.18Lg\n", i, wavefn[i], -logl(wavefn[i]));
	}
}

#endif

// Search for zero-crossings
void zero_cross(size_t len, double omega, double pow)
{
	long double *p = (long double*) malloc (len * sizeof(long double));
	for (size_t i=0; i<len; i++) p[i] = sigmaf(i, pow);
	set_pot (len, p);

	long double lo = 0.0;
	long double ylo = solve(len, omega, lo);
	for (long double hi=0.1; hi <30; hi+= 0.05)
	{
		long double yhi = solve(len, omega, hi);

		// If the two values are of opposite signes, then there's a sign crossing.
		if (ylo*yhi < 0.0)
		{
			long double eig = find_eig(len, omega, lo, hi);
			printf("pow=%g found %Lg\n", pow, eig);
		}
		lo = hi;
		ylo = yhi;
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

	double pow = atof(argv[2]);
	double omega = 1.0;

	init(len, omega);

	zero_cross(len, omega, pow);
}
