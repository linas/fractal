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

long double *wavefn;

long double *pot;

long double misma;

// Continuum limit for very small omega
long double omega = 0.08L;

void init(size_t len)
{
	size_t i;

	wavefn = (long double *) malloc(len * sizeof(long double));
	pot = (long double *) malloc(len * sizeof(long double));
	for (i=0; i<len; i++) wavefn[i] = 0.0L;

	// Simple harmonic oscillator!
	for (i=0; i<len; i++) pot[i] = 0.5L * i*i * omega*omega;

	// for (i=0; i<len; i++) pot[i] = totient_phi(i);
	// for (i=0; i<len; i++) pot[i] = divisor(i);
	// for (i=0; i<len; i++) pot[i] = sigma(i, 1);
	// for (i=0; i<len; i++) pot[i] = sigma(i, 2);
	// for (i=0; i<len; i++) pot[i] = sigma(i, 3);
	// for (i=0; i<len; i++) pot[i] = sigmalog(i, 2.0);
	// for (i=0; i<len; i++) pot[i] = sigmalog(i, 1.0);
	// for (i=0; i<len; i++) pot[i] = sigmaf(i, 1.2);

	// for (i=0; i<len; i++) printf("potential is %zd %Lf\n", i, pot[i]);

   wavefn[len-2] = 0.001;
}

// Solve teh schroedinger eqn 
// This resembles the bessel function, I guess ... 
void solve(size_t len, long double energy)
{
	energy *= omega;

	size_t i;
	for (i=len-2; 0 < i; i--)
	{
		long double ct = 2.0L * (1.0L + pot[i] - energy);
		wavefn[i-1] = ct * wavefn[i] - wavefn[i+1];
	}

	// Now renormalize
	long double acc = 0.0L;
	for (i=0; i<len; i++) acc += wavefn[i] * wavefn[i];

	acc = sqrtl(1.0L / acc);

	for (i=0; i<len; i++) wavefn[i] *= acc;

	misma = fabsl ((1.0L + pot[0] - energy) * wavefn[0] - 0.5 * wavefn[1]);
	misma -= fabsl (0.5 * wavefn[1]);
	// misma = fabsl(misma);
#if 0
	for (i=0; i<len; i++) 
		printf("wavefn= %zd %Lg\n", i, creall (wavefn[i]));
#endif
}

// Well, if it sbessel-like, then lets go with that idea...
void bessy(size_t len, double maxen)
{
	// printf("#\n#Bessel-like totient solution\n#\n");
	printf("#\n#Bessel-like divisor solution\n#\n");
	printf("# Energy  mismatch b0	b1	b2	b3	b4	b5\n");
	double en=0.0;
	double delta = maxen / 1000.0;
	for (en=0.0; en<maxen; en+=delta)
	{
		solve(len, en);

		printf("%g	%Lg	%Lg	%Lg	%Lg	%Lg	%Lg	%Lg\n",
			en, creall(misma),
			creall(wavefn[0]), creall(wavefn[1]),
			creall(wavefn[2]), creall(wavefn[3]),
			creall(wavefn[4]), creall(wavefn[5]));
	}
}

main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s <len> <maxenergy> <scale>\n", argv[0]);
		exit (1);
	}
	int len = atoi(argv[1]);
	printf("# Length=%d\n", len);

	double en = atof(argv[2]);
	omega = atof(argv[3]);

	init(len);

   // solve(len, en);
	bessy(len, en);
}
