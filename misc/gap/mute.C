
/* mute.C
 *
 * integrals of permutated, swapped continued fractions.
 *
 * Linas Vepstas October 2004
 */

#include <math.h>
#include <stdio.h>

#include "Farey.h"


double tegral (int p, int q, double t)
{
	ContinuedFraction f;
	// f.SetEvenize();

	int ngrid=12123;
	int i;
	double cacc = 0.0;
	double sacc = 0.0;
	for (i=0; i<ngrid; i++)
	{
		double x = ((double) (2*i+1))/((double) 2*ngrid);
		double term = cos (t*log (x))/sqrt(x);
		
		f.SetRatio (2*i+1,2*ngrid);
		f.SwapTerms (p,q);
		double y = f.ToReal();
		term *= y;
		cacc += term;
		
		term = sin (t*log (x))/sqrt(x);
		term *= y;
		sacc += term;
	}
	cacc /= ngrid;
	sacc /= ngrid;

	printf ("%8.6g	%8.6g	%8.6g\n", t, cacc, sacc);
	return cacc;
}

main (int argc, char *argv[])
{
	
	if (argc<3)
	{
		printf ("Usage: %s <p> <q>\n", argv[0]);
		exit (1);
	}
	int p = atoi (argv[1]);
	int q = atoi (argv[2]);
	printf ("# swapping terms %d and %d\n", p,q);

	double t=10.0;
	for (t=13.0; t<20; t +=0.01)
	{
		double y = tegral (p,q,t);


	}
}

