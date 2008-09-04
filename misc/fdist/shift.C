/*
 * shift.C
 *
 * Binary shift of continued fractions.
 *
 * Linas Vepstas Septemeber 2008
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

void prt_shift(int nbins, int kshift)
{
   ContinuedFraction f;
	for (int i=0; i<nbins; i++)
	{
		double x = ((double) i) / ((double) nbins);

   	f.SetRatio (i+1, nbins);
		f.BinaryLeftShift(kshift);
   	double y = f.ToReal (); 

#if 1
		printf ("%d	%8.6g	%8.6g\n", 
			i, x, y);
		fflush (stdout);
#endif
	}
}

main(int argc, char *argv[])
{
	int i;

	if (argc <3)
	{
		fprintf (stderr, "Usage: %s <nbins> <shift>\n", argv[0]);
		exit (1);
	}
	int nbins = atoi (argv[1]);
	int shift = atoi (argv[2]);

	prt_shift (nbins, shift);
}

