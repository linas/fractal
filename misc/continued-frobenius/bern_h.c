
/* bern_h.c
 *
 * Bernoulli hamiltonian
 *
 * Linas Vepstas November 2004
 */

#include "zetafn.h"
#include <stdio.h>
#include <math.h>

double elt (int j, int k)
{
	int n;
	double acc=0.0;

	for (n=j; n<= k; n++)
	{
		double term = binomial (n,j) * bernoulli (n-j);
		if (n != 0) term *= binomial (k,n-1);
		// term /= n;
		acc += term;
	}

	return acc;
}

int
main ()
{
	int j, k;

	for (j=0; j<7; j++) 
	{
		for (k=j; k<7; k++)
		{
			double e = elt (j,k);
			printf ("j=%d k=%d e=%g\n", j, k, e);
		}
		printf ("\n");
	}
	return 0;
}
