
/*
 * conj.c
 *
 * Conjugate map explorer
 * Goal: try to find conjugate maps of the tent map of a specific shape
 *
 * December 2004
 */

#include <math.h>
#include <stdlib.h>

double a = 0.5;

double 
phi (double x)
{
	double y;

	y = 0.5 + 0.5*a*(a+1)*(1-2*x)/(x*x-x-a*(a+1));

	return y;
}

double 
phiinv (double y)
{
	double rad;

	rad = 4*a*a*(a+1)*(a+1)/((2*y-1)*(2*y-1));
	rad += 1 + 4*a*(a+1);
	rad = 0.5*sqrt(rad);

	double b = 0.5 - a*(a+1)/(2*y-1);

	if (y<0.5) return b-rad;
	return b+rad;
}

double 
phit (double x)
{
	double y;
	double b;

	y = tan (a*M_PI*(x-0.5));
	b = 0.5/ tan (a*M_PI*0.5);
	y *= b;
	y +=0.5;

	return y;
}

double 
phitinv (double y)
{
	double b,x;

	b = 0.5/ tan (a*M_PI*0.5);
	x = atan ((y-0.5)/b);
	x /= a*M_PI;
	x +=0.5;

	return x;
}

double 
tent (double x)
{
	if (0.5 > x) return 2.0*x;
	return 2.0-2.0*x;
	// return 2.0*x-2.0;
}

double 
eff (double x)
{
	if (0.5 > x) return x/(1-x);
	return ((1-x)/x);
}
	
int
main(int argc, char ** argv)
{
	int i;

	a = atof (argv[1]);
	printf ("#  a=%g\n", a);

	int imax=233;
	for (i=0; i<imax; i++) 
	{
		double x = ((double) i)/((double) imax);

		double y = phitinv (x);
		y = tent (y);
		double xx = phit (y);

		printf ("%d %g  %g  %g	%g\n", i, x, y, xx, eff(x));
	}
	return 0;
}
