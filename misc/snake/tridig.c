
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

// Ternary (base-3) decomposition and reassembly
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

// Ternary (base-3) decomposition and middl-reflected reassembly
double reflect(double x, double w)
{

	double y = 0.0;
	double wn = w;

	for (int i=1; i<32; i++)
	{
		int dig = tridig(x); // ternary digit of x.

		if (1 == dig)
		{
			x = 2.0 - 3.0*x;     // right-shift
		}
		else
		{
			x = 3.0*x - dig;     // right-shift
		}
		y += dig * wn;
		wn *= w;
	}

	return y;
}

// Ternary (base-3) decomposition and zig-zag reassembly
double zigzag(double x, double w)
{

	double y = 0.0;
	double wn = w;

	for (int i=1; i<32; i++)
	{
		int dig = tridig(x); // ternary digit of x.
		x = 3.0*x - dig;     // right-shift

		if (0 == dig)
		{
			y += x * wn;
		}
		else if (1 == dig)
		{
			y += (1.0 - 2.0*x) * wn;
		}
		else
		{
			y += (x - 1.0) * wn;
		}
		wn *= w;
	}

	return y;
}

int main (int argc, char *argv[])
{
	int npts = 1800;
	double w = atof(argv[1]);
	for (int i=0; i<npts; i++)
	{
		double x = ((double) i) / ((double) npts);
		// double y = cantern(x, w);
		double y = zigzag(x, w);
		printf("%d	%g	%g\n", i,x,y);
	}
}
