//
// potential.cc
//
// Sanity check some numer-theoretic functions
//
// Linas Vepstas 23 Nov 2014
//
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "moebius.h"
#include "totient.h"

int divisor_sum(int n)
{
	int sum = 0;
	for (int i=0; i<=n; i++) sum += divisor(i);
	return sum;
}

main ()
{
	int len = 100;
	// For the simple harmonic oscillator, sqrt(omega) gives the length
	// scale.  len is the number of grid points.  We want the number of
	// grid points to extend deep into the classically forbidden region,
	// else the solver goes crazy.
	printf("# Length=%d\n", len);

	for (int i=0; i<len; i++)
	{
		// Euler totient function
		printf("%d	%g\n", i, (double) totient_phi(i));

		// Divisor function
		// printf("%d	%g\n", i, (double) divisor(i));
		// printf("%d	%g\n", i, (double) divisor_sum(i));

		// Sigma function (== divisor for power=0)
		// printf("%d	%g\n", i, (double) sigma(i, 1));
	}
}

