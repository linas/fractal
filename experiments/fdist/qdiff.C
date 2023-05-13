/*
 * qdiff.C
 *
 * Experimental diffs
 *
 * Linas March 2010
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"


main(int argc, char *argv[])
{
	int i;

	if (argc <2)
	{
		fprintf (stderr, "Usage: %s <epsilon>\n", argv[0]);
		exit (1);
	}
	double eps = atof(argv[1]);

   ContinuedFraction f;

	// double rt = 0.5 * sqrt(2.0);

	// Golden ratio
	double rt = 0.5 * (sqrt(5.0) - 1.0);
	rt = 1.0-rt;

  	// f.SetRatio (n, (unsigned long) RAND_MAX);
	f.SetReal(rt);
  	double cent = f.ToFarey (); 

	f.SetReal(rt-eps);
  	double lo = f.ToFarey (); 

	f.SetReal(rt+eps);
  	double hi = f.ToFarey (); 

	double dlo = (cent-lo) / eps;
	double dhi = (hi-cent) / eps;

	double disc = dhi-dlo;

	printf("its %g   %g   %g   disc=%g\n", cent, dlo, dhi, disc);
}

