
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
#include "totient.h"

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

void rebin (int nbins, int kmax)
{
	int i, k;

	double * rebin = (double*) malloc (nbins*sizeof(double));
	double * imbin = (double*) malloc (nbins*sizeof(double));

	// First compute the fourier xform.
	double resum = 0.0;
	double imsum = 0.0;
	for (i=0; i<nbins; i++)
	{
		/* Midpoint of the bin */
		double x = ((double) 2*i+1) / ((double) 2*nbins);

		double rev = 0.0;
		double imv = 0.0;
		for (k=0; k<kmax; k++)
		{
			double theta = 2.0*M_PI*x*k;
			double si = sin(theta);
			double co = cos(theta);

			// double fk = mangoldt_lambda (k);
			// double fk = totient_phi(k);
			// double fk = liouville_lambda (k);
			// double fk = liouville_omega (k);
			// double fk = mertens (k);
			double fk = moebius_mu(k);
			rev += co * fk;
			imv += si * fk;
		}

		resum += rev;
		imsum += imv;

		// printf ("%d	%g	%g	%g	%g	%g\n", i, x, rev, imv, resum, imsum);

		rebin[i] = rev;
		imbin[i] = imv;
	}

	// Next, undo the Fourier.
	for (k=0; k<kmax; k++)
	{
		double res = 0.0;
		double ims = 0.0;
		for (i=0; i<nbins; i++)
		{
			/* Midpoint of the bin */
			double x = ((double) 2*i+1) / ((double) 2*nbins);
			double theta = 2.0*M_PI*x*k;
			double si = sin(theta);
			double co = cos(theta);

			res += co * rebin[i];
			ims += si * rebin[i];
		}
		double sane = moebius_mu(k);
		printf ("%d	%g	%g	%g\n", k, res, ims, sane);
	}
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

	rebin (nbins, kmax);
}
