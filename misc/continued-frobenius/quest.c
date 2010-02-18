
/*
 * quest.c
 *
 * wild guess
 *
 * Linas Vepstas Feb 2010
 */

#include "question.h"

double gkw (double x)
{
	int n;
	double acc, term;

	acc = 0.0;
	for (n=1; n<500; n++)
	{
		term = 1.0 / (x + n);
		acc += term*term *question_inverse(1.0-term);
	}

	return acc;
}


main ()
{
	int i;
	double x, y;

	int nmax = 200;
	for (i=1; i<nmax; i++)
	{
		x = ((double) i) / ((double) nmax);
		y = gkw(x);

		printf ("%d	%f 	%f	%f\n", i, x, y, y/x);
	}

}
