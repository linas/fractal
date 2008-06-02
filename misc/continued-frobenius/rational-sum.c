/*
 * rational-sum.c
 *
 * rational zeta series sums.
 *
 * Linas June 2008
 */

#include <math.h>
#include <stdio.h>

#include "zetafn.h"

long double zlog2(void)
{
	int k;
	long double acc = 0.0L;
	for (k=1; k<60; k++)
	{
		long double term = zetam1 (2*k);
		term /= k;
		acc += term;
	}
	return acc;
}

int
main ()
{
	long double l2 = zlog2();
	printf ("its %Lg %Lg\n", l2, M_LN2-l2);

	return 0;
}
