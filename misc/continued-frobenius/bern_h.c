
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
		// term *= pow (0.5, n);
		acc += term;
	}

	return acc;
}

double basic_elt (int j, int k)
{
	if (j>k) return 0.0;
	if (j==k) return pow (0.5, k);
	return binomial(k,j)* pow (0.5, k+1);
}

int
main ()
{
	int j, k;

	int nmax = 23;
	for (j=1; j<nmax; j++) 
	{
		for (k=j; k<nmax; k++)
		{
			double e = elt (j,k);
			double b = basic_elt (j,k);
			// printf ("j=%d k=%d e=%g   \tb=%g\n", j, k, e, b);
			printf ("j=%d k=%d e=%15.10g  \n", j, k, e);
		}
		printf ("\n");
	}
	return 0;
}
