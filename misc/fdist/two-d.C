/*
 * Verify calculations of two-dimensional rep of dyadic group
 * and specifically verify eigenvalues of dyadic GKW
 *
 * March 2010
 */

#include <math.h>

double bee(double x)
{
	x -= floor(x);
	if (x<0.5) return 0.0;
	return 1.0;
}

double beew(double w, double x)
{
	int n;
	double wn = 1.0;
	double tn = 1.0;
	int nmax = floor(25.0 / log(w));
	double acc = 0.0;
	for (n=0; n<nmax; n++)
	{
		double term = wn * b(tn*x);
		acc +=term;
		wn *= w;
		tn *= 2.0;
	}
	return acc;
}
