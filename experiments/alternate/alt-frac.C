
/*
 * alt-frac.C
 *
 * Same idea as described in alt-binary, but applied to continued
 * fractions
 *
 * Linas Vepstas June 2014
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"


double falternate(int n, int d)
{
	ContinuedFraction f;
	f.SetRatio(n, d);

	double acc = 0.0;
	double sgn = 1.0;

	for (int i=1; i<f.GetNumTerms(); i++)
	{
		ContinuedFraction g = f;
		g.DropTerm(i);
		acc += sgn * g.ToReal();
		sgn = - sgn;
	}

	return acc;
}

main (int argc, char * argv[])
{
	int deno = 1513;
	for (int i=1; i<deno; i++)
	{
		double x = ((double) i) / ((double) deno);
		double y = falternate(i, deno);

		printf("%16.14g	%16.14g\n", x, y);
	}
}
