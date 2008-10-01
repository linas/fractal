
/*
 * Build a takaagi-style distribution
 *
 * Linas Vepstas October 2008
 */

#include <stdio.h>

double triangle (double x)
{
	x -= (int) x;
	if (x < 0.5) return x;
	return 1.0-x;
}

double prod (double x)
{
	int n;
	double prod = 1.0;
	double tk = 1.0;
	for (n=0; n< 20; n++)
	{
		prod *= 0.5 * (1.0 + triangle(tk*x));
		tk *=- 2.0;
	}
}

void graph(int npts)
{
	int i;

	double acc = 0.0;
	for (i=0; i<npts; i++)
	{
		double x = (double) i / ((double) npts);
		double y = prod(x);
		acc += y;

		printf ("%d	%g	%g %g\n", i, x, y, acc);
	}
}

main ()
{
	graph (567);
}
