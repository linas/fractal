
/*
 * tri.c
 *
 * ternaray decomposition
 * November 2017
 */

#include <stdio.h>

int tridig(double x)
{
	if (3.0*x < 1.0) return 0;
	if (3.0*x < 2.0) return 1;
	return 2;
}

int main ()
{

	double x = 4.0/7.0;
	double y = 0.0;
	double tn = 1.0/3.0;

	printf("its %g\n", x);
	for (int i=1; i<20; i++)
	{
		int dig = tridig(x);
		x = 3.0*x - dig;

		y += dig * tn;
		tn *= 1.0/3.0;

		printf("its %d %g %g %g\n", i,x,y,tn);
	}
}
