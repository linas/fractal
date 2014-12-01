//
// Quick double-check
//
// Linas Nov 2014

#include <math.h>
#include <stdio.h>

main()
{
	double x;
	for (x=0.0; x<1.0; x+=0.005)
	{
		double co=0.0;
		double si=0.0;
		for (int n=1; n<101001; n++)
		{
			co += cos(2.0*M_PI*n*x) / ((double) n);
			si -= sin(2.0*M_PI*n*x) / ((double) n);
		}

		co /= M_PI;
		si /= M_PI;

		// si should be Bernoulli polunomial B_1(x) = x-0.5
		// co should be Clausen function of order 1
		printf("%f	%f	%f	%f\n", x, co, si, x-0.5);
	}
}
