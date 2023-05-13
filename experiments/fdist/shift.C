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
#include "question.h"

void prt_shift(int nbins, int kshift)
{
   ContinuedFraction f;
	for (int i=0; i<nbins; i++)
	{
		double x = ((double) i) / ((double) nbins);

#if NOT_QUITE_IT
   	f.SetRatio (i+1, nbins);
		f.BinaryLeftShift(kshift);
   	double y = f.ToReal (); 
#endif
   	f.SetRatio (i+1, nbins);
		double far = f.ToFarey();
		for (int k=0; k<kshift; k++)
		{
			far *= 2.0;
		}
		far -= floor(far);
		if (0.5 < far) far = 1.0 - far;
		double y = question_inverse (far);

#if 1
		printf ("%d	%8.6g	%8.6g\n", 
			i, x, y);
		fflush (stdout);
#endif
	}
}

void prt_qprime(int nbins)
{
   ContinuedFraction f;

	double gral = 0.0;
	for (int i=0; i<nbins; i++)
	{
		double x = ((double) i) / ((double) nbins);

   	f.SetRatio (i+1, nbins);
		double qprime = 1.0;
		for (int kshift = 0; kshift < 20; kshift++)
		{
			f.BinaryLeftShift(kshift);
   		double y = f.ToReal (); 
			qprime /= (1.0-y)*(1.0-y);
		}
		gral += qprime;

#if 1
		printf ("%d	%8.6g	%8.6g	%8.6g\n", 
			i, x, qprime, gral);
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
	// prt_qprime (nbins);
}

