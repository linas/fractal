
/*
 * similarity -conjugacy
 *
 * Linas Vepstas Feb 2010
 */

#include <math.h>
#include <stdio.h>

#include "question.h"

double ach(double x)
{
	double y = 1.0/x;
	y -= floor(y);
	return y;
}

double cee(double x)
{
	int tk;
	double y = 1.0/x;
	int tn = (int) floor(y);
	if (tn < 0) tn  = 0;

	tk = 1;
	while(tn !=0)
	{
		tn >>= 1;
		tk <<= 1;
	}
	y = 2.0 - tk*x;
	return y;
}

double funky(double x)
{
	double y;
	y = question_inverse(x);
	// y = 2.0*y-1.0;
	y = 2.0*y-1.0;
	y = fquestion_mark(y);
	return y;
}

double gkw (double x)
{
	int k;
	double term, acc;

	acc = 0.0;
	for (k=1; k<500; k++)
	{
		term = 1.0 / (x+ ((double)k));
		acc += term*term*funky(term);
	}

	return acc;
}

int main (int argc, char * argv[])
{
	int npts = 300;
	int i;

	npts = atoi(argv[1]);

	for (i=0; i<npts; i++)
	{
		double x = ((double) i) / ((double) npts);
#if 0
		double g;
		double y = funky(x);
		g = gkw(x);
		printf ("%d	%f	%f	%f\n", i, x, y, g);
#endif

		double h, c;
		h = question_inverse(x);
		h = ach(h);
		h = fquestion_mark(h);

		c = x;
		c = cee(c);

		printf("%d	%f	%f	%f	%g\n", i, x, c, h, c-h);
	}

	return 0;
}
