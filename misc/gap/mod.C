
/* 
 * explore limits of modular functions 
 *
 * Linas February 2005 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "modular.h"

double j_limit (double x, double r)
{
	double jre, jim;

	double qre, qim;
	qre = r*cos (2.0*M_PI*x);
	qim = r*sin (2.0*M_PI*x);
	klein_j_invariant_c (qre, qim, &jre, &jim);

	double jabs = sqrt (jre*jre +jim*jim);
	jabs = 1.0/jabs;
	return jabs;
}

main (int argc, char * argv[]) 
{
	int nmax = 400;
	int i;

	double r = 0.3;

	r = atof (argv[1]);

	double acc = 0.0;
	for (i=0; i<nmax; i++)
	{
		double x = ((double) i)/ ((double) nmax);

		double y = j_limit (x, r);
		acc += y;
		printf ("%d	%g	%g	%g\n", i, x, y, acc);
	}
}
