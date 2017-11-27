
/*
 * tri.c
 *
 * ternaray decomposition
 * November 2017
 */

#include <stdio.h>
#include <stdlib.h>

int tridig(double x)
{
	if (3.0*x < 1.0) return 0;
	if (3.0*x < 2.0) return 1;
	return 2;
}

double cantern(double x, double w)
{

	double y = 0.0;
	double wn = w;

	for (int i=1; i<32; i++)
	{
		int dig = tridig(x); // ternary digit of x.
		x = 3.0*x - dig;     // right-shift

		y += dig * wn;
		wn *= w;
	}

	return y;
}

int main (int argc, char *argv[])
{
	int npts = 400;
	double w = atof(argv[1]);
	for (int i=0; i<npts; i++)
	{
		double x = ((double) i) / ((double) npts);
		double y = cantern(x, w);
		printf("%d	%g	%g\n", i,x,y);
	}
}
