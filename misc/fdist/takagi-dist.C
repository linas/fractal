
/*
 * Build a takaagi-style distribution
 *
 * Linas Vepstas October 2008
 */

#include <math.h>
#include <stdio.h>

double triangle (double x)
{
	x -= (int) x;
	if (x < 0.5) return x;
	return 1.0-x;
}

double a_k (double x, int k)
{
	int i;
	for (i=0; i<k; i++)
	{
		if (x < 0.5) x = x/(1.0-x);
		else x = (2.0*x-1)/x; 
		// else x = (1.0-x)/x; 
	}
	return triangle (x);
}

double prod (double x)
{
	int n;
	double prod = 1.0;
	double tk = 1.0;
	for (n=0; n< 20; n++)
	{
		// tk *= 2.0;
		// double mand = 1.0 + triangle(tk*x);
		// double mand =  1.25 + 0.25*sin(M_PI*tk*x);
		double mand = 1.0 + a_k(x,n+1);
		// prod *= 0.5 * mand*mand;
		prod /= 0.5 * mand*mand;
	}

	return prod;
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
	graph (1567);
}
