//
// Quikc double-check
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
		printf("%f	%f	%f	%f\n", x, co, si, x-0.5);
	}
}
