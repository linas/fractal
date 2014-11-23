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

#include "moebius.h"
#include "totient.h"

#define LEN 130
long double complex wavefn[LEN];

long double pot[LEN];

void init()
{
	int i;
	for (i=0; i<LEN; i++) wavefn[i] = 0.0;
	for (i=0; i<LEN; i++) pot[i] = totient_phi(i);
	for (i=0; i<LEN; i++) printf("potential is %d %Lf\n", i, pot[i]);

   wavefn[LEN-2] = 0.001;
}

void solve(long double energy)
{
	int i;
	for (i=LEN-2; 0 < i; i--)
	{
		long double ct = 2.0 + pot[i] - energy;
		wavefn[i-1] = ct * wavefn[i] - wavefn[i+1];
	printf("duude fn= %Lg\n", creall (wavefn[i-1]));
	}
}

main()
{
	init();
   solve(9.0);
}
