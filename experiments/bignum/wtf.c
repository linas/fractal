
#include <math.h>
#include <stdio.h>

double summy (int k)
{
	double acc = 0.0;
	int n;

	for (n=1; n<40; n++)
	{
		double term  = exp (2.0*M_PI*n) - 1.0;
		term *= n*n*n;
		acc += 1.0/term;
	}

	return acc;
}

main ()
{
	double x = summy(3);
	printf ("its %g\n", x);
}
