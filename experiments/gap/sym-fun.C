
/* sym-fun.C
 *
 * Sawtooth and related functions generated by the Minkowski symmetery group.
 *
 * Linas Vepstas October 2004
 */

#include <math.h>
#include <stdio.h>

#include "Farey.h"
#include "question.h"

main (int argc, char *argv[])
{
	int i;
	ContinuedFraction f;
	f.SetEvenize();


	int nmax = 1431;
	for (i=0; i<nmax; i++)
	{
		double x = ((double) 2*i+1)/ ((double) (2*nmax));

#define SAWTOOTH
#ifdef SAWTOOTH
		// this is what continued-fraction truncation looks like 
		// in binary-digit space
		double z = question_inverse(x);
		double y = 1.0/z;
		y -= floor (y);
		f.SetReal (y);
		double qy = f.ToFarey();
		y = qy;
#endif

#ifdef BERNOULLI
		// this is what binary-digit truncation looks like in
		// in continued fraction space.
		f.SetReal (x);
		double yb = f.ToFarey();
		yb *= 2.0;
		yb -= floor (yb);
		double y = question_inverse (yb);
#endif

#if EIGENVALUE_ZERO
		// This is the zeroth eigenvalue, in binary-space
		double xc = question_inverse(x);
		double yc = 1.0/(1.0+xc);
		f.SetReal (yc);
		double y = f.ToFarey();
		double gy = 1.0 - 0.5*x;
#endif

#ifdef WHATEVER
		// This is the 1st approx eigenvalue, in binary-space
		double xc = question_inverse(x);
		// double yc = 1.0 / ((1.0+xc)*(1.0+xc)*sqrt(1.0+xc));
		double yc = (2.0*xc-1.0)/(2.0*xc+1.0);
		f.SetReal (yc);
		double y = f.ToFarey();
#endif

		double gy = 1.0 - 0.5*x;
		printf ("%5d	%8.6g	%8.6g	%8.6g\n", i, x, y, y-gy);

	}
}

