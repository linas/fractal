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

/// The continuum limit is obtained for very small omega.
/// Note that, for the simple harmonic oscillator, the length scale
/// is given by sqrt(mass*omega). Here we take mass=1.0

void init(size_t len, long double omega)
{
	// For the simple harmonic oscillator, sqrt(omega) gives the length
	// scale.  len is the number of grid points.  We want the number of
	// grid points to extend deep into the classically forbidden region,
	// else the solver goes crazy.
	len /= sqrt(omega);
	printf("# Length=%zd omega=%Lg\n", len, omega);

	wavefn = (long double *) malloc(len * sizeof(long double));
	pot = (long double *) malloc(len * sizeof(long double));

	size_t i;

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
}

/// Solve teh schroedinger eqn.
///
/// 'len' is the number of grid points. It should be large enough so
/// that the wave-function starts out deep in the classically-forbidden
/// zone (i.e. it is assumed that the potential energy is mostly
/// positive, and increasing, so that in most spots far from the origin
/// is is much larger than the energy.)  The solution technique here
/// is thus very similar to that for spherical bessel functions: one
/// starts deep in the forbidden region, and backwards-iterates. This
/// causes the error to be suppressed geometrically with each iteration.
///
/// The length scale is given by sqrt(mass*omega), which is exactly
/// correct for the simple harmonic oscillator.  We take mass=1 here.
/// The corresponding energy scale is simply 1/omega.
/// 
/// The argument 'energy' is the guess for the wave-function eigenvalue.
/// (It is rescaled by omega to get the grid units).
///
/// The returned value is the discontinuity at the origin. This can be
/// used as a boundary condition: the discontinuity must be zero to
/// get true eigenstates of the double-sided symmetric oscillator.
///
/// Another possible boundary condition is to set the wave-function to
/// zero at the origin (i.e. a hard-wall boundary condition).
///
long double solve(size_t len, long double omega, long double energy)
{
	// Rescale to work in grid units.  This is transparent to the user,
	// who specifies only the 'physical' energy.
	energy *= omega;

	// Laplacian is -psi(i-1) + 2*psi(i) - psi(i+1)
	// Schreodinger eqn is 1/2 Lap(psi) + V(i) psi = E psi
   wavefn[len-1] = 0.0L;
   wavefn[len-2] = 1.0e-100L;
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

	long double discon = fabsl ((1.0L + pot[0] - energy) * wavefn[0] - 0.5 * wavefn[1]);
	discon -= fabsl (0.5 * wavefn[1]);
	// discon = fabsl(discon);
#if 0
	for (i=0; i<len; i++) 
		printf("wavefn= %zd %Lg\n", i, creall (wavefn[i]));
#endif

	return discon;
}

// Well, if its bessel-like, then lets go with that idea...
void bessy(size_t len, double omega, double maxen)
{
	// printf("#\n#Bessel-like totient solution\n#\n");
	printf("#\n# Divisor solution\n#\n");
	printf("# Energy  discontinuity b0	b1	b2	b3	b4	b5\n");
	double en=0.0;
	double delta = maxen / 1000.0;
	for (en=0.0; en<maxen; en+=delta)
	{
		long double discon = solve(len, omega, en);

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
