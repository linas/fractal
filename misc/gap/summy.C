
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

main () 
{
	int i;

	int nmax = 23;

	for (i=0; i<nmax; i++)
	{
		double x = ((double) i ) / ((double) nmax);
		double y = summy (x);

		double z = 1.0;

		printf ("%8.6g	%8.6g	%8.6g\n", x,y,z);
	}

}
