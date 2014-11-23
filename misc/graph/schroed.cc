//
// schroed.cc
//
// Solve schroedinger eqn ion a 1D lattice.
//
// Linas Vepstas 23 Nov 2014
//
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "moebius.h"
#include "totient.h"

long double complex *wavefn;

long double *pot;

void init(size_t len)
{
	size_t i;

	wavefn = (long double complex *) malloc(len * sizeof(long double complex));
	pot = (long double *) malloc(len * sizeof(long double));
	for (i=0; i<len; i++) wavefn[i] = 0.0;
	for (i=0; i<len; i++) pot[i] = totient_phi(i);
	for (i=0; i<len; i++) printf("potential is %zd %Lf\n", i, pot[i]);

   wavefn[len-2] = 0.001;
}

void solve(size_t len, long double energy)
{
	size_t i;
	for (i=len-2; 0 < i; i--)
	{
		long double ct = 2.0 + pot[i] - energy;
		wavefn[i-1] = ct * wavefn[i] - wavefn[i+1];
	}

	// Now renormalize
	long double acc = 0.0;
	for (i=0; i<len; i++) acc += creall(wavefn[i] * conjl(wavefn[i]));

	acc = sqrtl(1.0 / acc);

	for (i=0; i<len; i++) wavefn[i] *= acc;
	for (i=0; i<len; i++) 
		printf("wavefn= %zd %Lg\n", i, creall (wavefn[i]));
}

main(int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <len>\n", argv[0]);
		exit (1);
	}
	int len = atoi(argv[1]);
	printf("Length=%d\n", len);

	init(len);
   solve(len, 9.0);
}
