
/*
 * quest.c
 *
 * wild guess
 *
 * Linas Vepstas Feb 2010
 */

#include <math.h>
#include <stdio.h>

#include "question.h"

double gkw (double x)
{
	int n;
	double acc, term;

	acc = 0.0;
	for (n=1; n<2500; n++)
	{
		term = 1.0 / (x + n);
		acc += term*term *question_inverse(1.0-term);
	}

	return acc;
}


int
main ()
{
	int i;
	double x, y, w;

	int nmax = 200;
	for (i=1; i<nmax; i++)
	{
		x = ((double) i) / ((double) nmax);
		w = question_inverse(1.0-x);
		y = gkw(x);

		printf ("%d	%f %f	%f	%f\n", i, x, w, y, y/w);
	}
	return 0;
}
