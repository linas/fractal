//
// schroed.cc
//
// Solve schroedinger eqn in a 1D lattice.  Focus on number-theoretic
// potentials.
//
// Linas Vepstas 23 Nov 2014
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "moebius.h"
#include "totient.h"

#include "schroed.h"

long double *wavefn = NULL;
long double *pot = NULL;

int enscale = 1;

int divisor_sum(int n)
{
	int sum = 0;
	for (int i=0; i<=n; i++) sum += divisor(i);
	return sum;
}

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
	// enscale = 1;  // for the sho
	// for (i=0; i<len; i++) pot[i] = 0.5L * i*i * omega*omega;

	enscale = 0;  // disable -- the number-theortic functions are scale-free

	// Euler totient function
	for (i=0; i<len; i++) pot[i] = totient_phi(i);

	// Divisor function
	// for (i=0; i<len; i++) pot[i] = divisor(i);
	// for (i=0; i<len; i++) pot[i] = divisor_sum(i);

	// Sigma function (== divisor for power=0)
	// for (i=0; i<len; i++) pot[i] = sigma(i, 1);
	// for (i=0; i<len; i++) pot[i] = sigma(i, 2);
	// for (i=0; i<len; i++) pot[i] = sigma(i, 3);
	// for (i=0; i<len; i++) pot[i] = sigmalog(i, 2.0);
	// for (i=0; i<len; i++) pot[i] = sigmalog(i, 1.0);
	// for (i=0; i<len; i++) pot[i] = sigmaf(i, 1.2);
	// for (i=0; i<len; i++) pot[i] = sigmaf(i, 2.3456);

	// for (i=0; i<len; i++) printf("potential is %zd %Lf\n", i, pot[i]);
}

void set_pot(size_t len, long double *p)
{
	if (pot) free(pot);
	pot = p;
}

/// Solve the schroedinger eqn.
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
	for (int j = 0; j<enscale; j++) energy *= omega;

	// Laplacian is -psi(i-1) + 2*psi(i) - psi(i+1)
	// Schroedinger eqn is 1/2 Lap(psi) + V(i) psi = E psi
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

#define REFLECT_AT_ZERO 1
#if REFLECT_AT_ZERO
	long double discon = fabsl ((1.0L + pot[0] - energy) * wavefn[0] - 0.5 * wavefn[1]);
	discon -= fabsl (0.5 * wavefn[1]);
#else

	long double discon = fabsl ((1.0L + pot[1] - energy) * wavefn[1] - 0.5 * wavefn[2]);
	discon -= fabsl (0.5 * wavefn[2]);
#endif
	// discon = fabsl(discon);
#if 0
	for (i=0; i<len; i++) 
		printf("wavefn= %zd %Lg\n", i, creall (wavefn[i]));
#endif

	return discon;
}

// ===============================================================

// Find zero. (aka find the eigenvalue). The best-guess eigenvalue
// should be brakceted by hi, lo.
//
// If 'odd' is set, then we look for eigenfunctions that are odd
// at the origin (i.e. the eigenfunction is zero at the origin)
// If it is not set, then we are looking for eignefunctions that
// are even (i.e. that the slope of the wave-function is zero at
// the origin)  Currently, solve returns that slope .. err, actually,
// it returns twice the slope; the 'discontinuity' of the second
// derivative.
long double find_eig(size_t len, long double omega, bool odd, long double lo, long double hi)
{
	long double dilo = solve(len, omega, lo);
	if (odd) dilo = wavefn[0];
	long double dihi = solve(len, omega, hi);
	if (odd) dihi = wavefn[0];

	int iter = 0;
	long double xmid, ymid;
	for (iter=0; iter<60; iter++)
	{
		long double delta = hi - lo;
		if (fabsl(delta) < 1.0e-17) return 0.5L*(hi+lo);

		// compute the midpoint
		long double slope = (dihi - dilo) / (hi - lo);

		xmid = lo - (dilo / slope);
		ymid = solve(len, omega, xmid);
		if (odd) ymid = wavefn[0];

#if 1
		printf("# duuude lo %Lg %Lg\n", lo, dilo);
		printf("# duuude mid %Lg %Lg\n", xmid, ymid);
		printf("# duuude hi %Lg %Lg\n#\n", hi, dihi);
#endif

		// secant method -- shuffle down by one
		lo = hi;
		dilo = dihi;
		hi = xmid;
		dihi = ymid;
	}

	printf("Error: failed to converge, ymid=%Lg at %21.18Lg\n", ymid, xmid);
	exit(1);
}

