
/*
 * similarity -conjugacy
 *
 * Linas Vepstas Feb 2010
 */

#include <math.h>
#include <stdio.h>

#include "question.h"

double funky(double x)
{
	double y;
	y = question_inverse(x);
	y = 2.0*y-1.0;
	y = fquestion_mark(y);
	return y;
}

double gkw (double x)
{
	int k;
	double term, acc;

	acc = 0.0;
	for (k=1; k<1500; k++)
	{
		term = 1.0 / (x+k);
		acc += term*term*funky(term);
	}

	return acc;
}

main ()
{
	int npts = 300;
	int i;

	for (i=0; i<npts; i++)
	{
		double g;
		double x = ((double) i) / ((double) npts);
		double y = funky(x);
		g = gkw(x);

		printf ("%d	%f	%f	%f\n", i, x, y, g);
	}

	return 0;
}
