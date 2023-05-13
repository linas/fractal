
/* integ.c
 * 
 * quick cheap cheesy integrator
 */

#include <math.h>

inline double 
m0 (double x)
{
	double v = 1.0/(1.0 +x);
	return (v * sqrt(8.0 / 3.0));
}

inline double 
m1 (double x)
{
	double y = 1.0 / (1.0+x);
	double v = y*y - 7.0 *y / 9.0;
	v *= 72.0 / sqrt (39.0);
	return v;
}

inline double 
m2 (double x)
{
	double y = 1.0 / (1.0+x);
	double v = y*y*y - 99.0*y*y/65.0 + 291.0*y/(8.0*65.0);
	v *= 40.0* sqrt ( 8.0 * 13.0 / 21.0);
	return v;
}

main () 
{
	double a00 = 0.0;
	double a01 = 0.0;
	double a02 = 0.0;
	double a11 = 0.0;
	double a12 = 0.0;
	double a22 = 0.0;
#define NMAX 12110123

	int i;
	for (i=0; i<NMAX; i++)
	{
		double x = ((double) i) / ((double) NMAX);
		double w = 1.0 / (1.0+x);

		a00 += w * m0(x) * m0(x);
		a01 += w * m0(x) * m1(x);
		a02 += w * m0(x) * m2(x);
		a11 += w * m1(x) * m1(x);
		a12 += w * m1(x) * m2(x);
		a22 += w * m2(x) * m2(x);
	}
	a00 /= ((double) NMAX);
	a01 /= ((double) NMAX);
	a02 /= ((double) NMAX);
	a11 /= ((double) NMAX);
	a12 /= ((double) NMAX);
	a22 /= ((double) NMAX);

	printf ("results a00=%f a11=%f a22=%f \n", a00, a11, a22);
	printf ("results a01=%f a02=%f a12=%f \n", a01, a02, a12);
}
