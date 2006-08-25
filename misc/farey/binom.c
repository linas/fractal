
/*
 * binom.c
 *
 * Experiments with binomial coefficients
 *
 * August 2006 Linas
 */

#include <stdio.h>
#include "binomial.h"

main()
{
	int n=30;

	int k;
	for (k=1; k<=n; k++)
	{
		int g = gcf32 (n,k);
		int p = n/g;
		int q = k/g;
		double b = fbinomial (p,q);
		double x = ((double)q)/((double)p);
		printf ("%g	%d	%d	%g\n", x, q,p,b);
	}
	
}
