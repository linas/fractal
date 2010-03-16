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

	if (argc <1)
	{
		fprintf (stderr, "Usage: %s\n", argv[0]);
		exit (1);
	}

   ContinuedFraction f;

	double rt = 0.5 * sqrt(2.0);
  	// f.SetRatio (n, (unsigned long) RAND_MAX);
	f.set(rt);
  	double x = f.ToFarey (); 

	printf("its %g\n", x);
}

