/*
 * fanal.C
 * 
 * Fourier transform of
 * Distribution of the Farey Numbers on the unit interval
 *
 * Linas October 2004
 * Linas Sept 2008
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

void fourier (int npts, int freq_max)
{
	int i, n;

	double gral = 0.0;

	/* Compute the integral of the distribution */
   ContinuedFraction f;

	for(i=0; i<npts; i++)
	{
		double x = ((double) i) / ((double) npts);

		/* Likewise, the midpoint */
   	f.SetRatio (2*i+1, 2*npts);
   	double far = f.ToFarey (); 

		double fx = 18.0;
		double tp = 1.0;
		for (n=1; n<freq_max; n++)
		{
			fx -= 2.0*cos(2.0*M_PI*tp*far);
			tp *= 2.0;
		}
		gral += fx / ((double) npts);
		printf ("%8.6g	%8.6g	%8.6g	%8.6g\n", x, fx, gral, far);
	}

}

main(int argc, char *argv[])
{
	int i;

	if (argc < 3)
	{
		fprintf (stderr, "Usage: %s <nbins> <freq>\n", argv[0]);
		exit (1);
	}
	int nbins = atoi (argv[1]);
	int max_freq = atoi (argv[2]);

	fourier (nbins, max_freq);
}

