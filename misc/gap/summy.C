
/*
 * summy.C
 * whacky sums of the GKW operator
 *
 * Linas Vepstas October 2004
 */

#include <math.h>

#include "Farey.h"

double 
summy (double x)
{
	ContinuedFraction f;
	f.SetEvenize();

	f.SetReal (x);

	double can = 2.0 - f.ToFarey();

	int n;
	double tn = 2.0;
	double acc = 0.0;
	for (n=1; n<20; n++)
	{
		double term = (can + tn) * (n+x)*(n+x);
		acc += 1.0/term;
		tn *= 2.0;
	}
	return acc*can;
}

double 
sum4 (double x)
{

	int n;
	double tn = 2.0;
	double acc = 0.0;
	for (n=1; n<20; n++)
	{
		double term = 2.0/(x+n) - 1.0 / (x+n+1);
		term *= term;
		term /= tn;
		acc += term;
		tn *= 2.0;
	}

	acc -= 8.0;
	acc *= 0.25;
	return acc;
}


double 
sum1 (double x)
{

	int n;
	double tn = 2.0;
	double acc = 0.0;
	for (n=1; n<20; n++)
	{
		double term = 1.0/(x+n);
		term -= 1.0/(x+n+1);
		term /=  tn;
		acc += term;
		tn *= 2.0;
	}

	return acc;
}

double 
sum2 (double x)
{

	int n;
	double tn = 4.0;
	double acc = 0.0;
	for (n=1; n<20; n++)
	{
		double term = 1.0/((x+n)*(x+n));
		term /= tn;
		acc += term;
		tn *= 4.0;
	}

	acc = (8.0*acc -1.0);
	return acc;
}

double gkw ( double (*fn)(double), double x)
{
	int n;

	double acc = 0.0;
	for (n=1; n<5000; n++)
	{
		double y = 1.0 / (x+((double)n));
		acc += fn(y) *y*y;
	}
	return acc;
}

main () 
{
	int i;

	int nmax = 523;

	double gm = 0.5* (sqrt(5) -1.0);

	for (i=0; i<nmax; i++)
	{
		double x = ((double) i ) / ((double) nmax);
		// double y = summy (x);

		double y = sum1(x);
		double z = 1.0;
		// double z = gkw (sum2, x);

		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x,y,z);
	}

}
