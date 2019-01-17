/*
 * gpf-dirichlet.c
 *
 * Misc assorted dircihlet sums and Mobius inversions of GPF.
 * January 2019
 */

#include <math.h>
#include <stdio.h>

#include "moebius.h"

/// Sum over divisors.
/// Same as Dirichlet convolution with the unit function.
long divisor_sum (long n, long (*f)(long))
{
	long sum = 0;
	for (long d=1; 2*d <= n; d++)
	{
		if (n%d) continue;
		sum += f(d);
	}
	sum += f(n);
	return sum;
}

/// Mobius inversion
/// Same as Dirichlet convolution with the Mobius function.
long moebius_invert (long n, long (*f)(long))
{
	long sum = 0;
	for (long d=1; 2*d <= n; d++)
	{
		if (n%d) continue;
		sum += f(d) * moebius_mu (n/d);
	}
	sum += f(n);
	return sum;
}

long unit (long n) { return 1; }

int main(int argc, char* argv[])
{
	for (long n=1; n<25; n++)
	{
		long sum = divisor_sum(n, unit);
		long inv = moebius_invert(n, unit);
		long gfn = unit(n);
		printf("n=%ld	gpf=%ld	sum=%ld	inv=%ld\n", n, gfn, sum, inv);
	}
}
