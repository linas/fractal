
/*
 * fiddle.C
 *
 * misc fiddleing
 */

#include <math.h>
#include <stdio.h>

#include "Farey.h"


double 
saw_eigenvalue (int m)
{
	int i, n;
	long double acc = 0.0;
	for (n=1;n<1000000; n++)
	{
		long double nnp = n*(n+1);
		long double term = n*(n+1);
		for (i=0; i<m; i++)
		{
			term *= nnp;
		}
		term = 1.0 / term;
		if (1.0e-16> term) break;
		acc += term;
	}

	return acc;
}

main ()
{
	int i;
	ContinuedFraction f;

	double z2 = M_PI*M_PI/6.0;
	double l1 = 2.0-z2;

	f.SetReal (l1);
	double c1 = f.ToFarey();

	printf ("its ?(%g)=%g\n", l1, c1);

	double acc = 0.0;
	double cacc = 0.0;
	double tn = 2.0;
	for (i=1; i<20; i++)
	{
		double eig = saw_eigenvalue (i);
		acc += eig;
		f.SetReal (eig);
		double qeig = f.ToFarey();
		cacc += qeig;
		printf ("%d    %g (x2=%g)  acc=%g   ?(eig)=%g   sum=%g\n", i, eig, tn*eig, acc, qeig, cacc);
		tn *= 2.0;
	}
	
}

