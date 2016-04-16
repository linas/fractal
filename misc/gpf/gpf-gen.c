/*
 * Generating functions for greatest prime factors.
 *
 * April 2016
 */

#include <math.h>
#include <stdio.h>

#include <gpf.h>

/*
 * Ordinary generating function for the greatest common factor.
 */
double gpf_ordinary(double x)
{
	double sum = 0;
	double xn = x;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn;
		xn *= x;
		if (n*xn < 1.0e-16) break;
	}

	return sum;
}

int main()
{
	for (double x=0.0; x< 1.0; x+= 0.002)
	{
		double y = gpf_ordinary(x);
		printf("%g\t%g\n", x, y);
	}
}
