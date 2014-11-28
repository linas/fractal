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
#include "totient.h"

#include "schroed.h"

// Search for zero-crossings.  The eigenvalue solving the schroed eqn
// is the one where the 'solve' function returns zero.
void zero_cross(size_t len, double omega, double pow)
{
	long double *p = (long double*) malloc (len * sizeof(long double));
	// for (size_t i=0; i<len; i++) p[i] = sigmaf(i, pow);
	for (size_t i=0; i<len; i++) p[i] = totient_phi(i);
	set_pot (len, p);

	bool odd = false;
	int nfound = 0;

	long double lo = 0.0;
	long double ylo = solve(len, omega, lo);
	for (long double hi=0.1; hi <30; hi+= 0.05)
	{
		long double yhi = solve(len, omega, hi);
		printf("pow=%g hi=%Lg yhi=%Lg\n", pow, hi, yhi);

		// If the two values are of opposite signes, then there's a sign crossing.
		if (ylo*yhi < 0.0)
		{
			long double eig = find_eig(len, omega, odd, lo, hi);
			nfound ++;
			printf("n=%d pow=%g found %20.17Lg\n", nfound, pow, eig);
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
