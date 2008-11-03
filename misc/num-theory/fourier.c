
/*
 * fourier.c
 *
 * Fourier transforms over varius number theoretic funcions
 *
 * Linas Vepstas November 2008
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "moebius.h"

int mertens (int n)
{
	int acc = 0.0;

	int i=1;
	for (i=1; i<=n; i++)
	{
		acc += moebius_mu (i);
	}

	return acc;
}

int main (int argc, char * argv[])
{
	int i, k;

	if (argc < 3)
	{
		fprintf (stdout, "Usage: %s <nbins> <fmax>\n", argv[0]);
		exit(1);
	}
	int nbins = atoi (argv[1]);
	int kmax = atoi(argv[2]);

	for (i=0; i<nbins; i++)
	{
		double x = ((double) 2*i+1) / ((double) 2*nbins);

		double reacc = 0.0;
		double imacc = 0.0;
		for (k=0; k<kmax; k++)
		{
			double si = sin(2.0*M_PI*x*k);
			double co = cos(2.0*M_PI*x*k);

			// mert = mertens (i);
			reacc += co * moebius_mu(i);
			imacc += si * moebius_mu(i);
		}

		printf ("%d	%g	%g	%g\n", i, x, reacc, imacc);
	}
}
