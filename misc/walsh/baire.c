/*
 * Rademacher functions for Baire space
 *
 * Used to illustrate the gkw.lyx paper.
 * Linas January 2016
 */

#include <math.h>
#include <stdio.h>

double h(double x)
{
	double y = 1.0 / x;
	return y - floor(y);
}

double rade(int n, int b, double x)
{
	for (int i=1; i<n; i++)
	{
		x = h(x);
	}
	x = 1.0 / x;
	if (b < x && x <= (b+1.0)) return 1.0;
	return 0.0;
}

int main(int argc, char* argv[])
{
	double x;
	for (x=0.0; x<1.0; x+=0.0001)
	{
		double r1 = rade(1,3,x);
		double r2 = rade(2,1,x);
		double r3 = rade(3,1,x);
		printf("%f	%f	%f	%f\n", x, r1, r2, r3);
	}

	return 0;
}

