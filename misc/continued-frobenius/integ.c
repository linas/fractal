
/* integ.c
 * 
 * quick cheap cheesy integrator
 */

#include <math.h>

double 
m0 (double x)
{
	double v = 1.0/(1.0 +x);
	return (v * sqrt(8.0 / 3.0));
}

double 
m1 (double x)
{
	double v = - 7.0 / (9.0 * (1.0 +x));
	v += 1.0 / ((1.0+x) *(1.0+x));
	v *= 72.0 / sqrt (39.0);
	return v;
}

main () 
{
	double a00 = 0.0;
	double a01 = 0.0;
	double a11 = 0.0;
#define NMAX 1100123

	int i;
	for (i=0; i<NMAX; i++)
	{
		double x = ((double) i) / ((double) NMAX);
		double w = 1.0 / (1.0+x);

		a00 += w * m0(x) * m0(x);
		a01 += w * m0(x) * m1(x);
		a11 += w * m1(x) * m1(x);
	}
	a00 /= ((double) NMAX);
	a01 /= ((double) NMAX);
	a11 /= ((double) NMAX);

	printf ("results %f %f %f \n", a00, a01, a11);
}
