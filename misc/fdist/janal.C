/*
 * janal.C
 * 
 * Fourier transform of
 * Distribution of the Farey Numbers on the unit interval
 * From hypothesized first principles.
 *
 * Linas October 2004
 * Linas Sept 2008
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double *bin = NULL;
int nbins = 0;

void alloc_bins(int nb)
{
	nbins = nb;
	bin = (double *) malloc(nb * sizeof(double));
	for (int i=0; i<nb; i++)
	{
		bin[i] = 0.0;
	}
}

void prt_bins(void)
{
	for (int i=0; i<nbins; i++)
	{
		double x = ((double) i) / ((double) nbins);
		printf ("%8.6g	%8.6g\n", x, bin[i]);
	}
}

void add_to_bin(double position, double amount)
{
	// First subtact int value of x
	// position -= floor (position);

	// position*nbins then truncate and bincount.
	double x = position * nbins;
	int i = (int) x;
	bin[i] += amount;

#ifdef MIRROR
	// mirror image
	position = 1.0 - position;
	x = position * nbins;
	i = (int) x;
	bin[i] += amount;
#endif
}

void hseq (double rlo, double rhi, double weight, double scale)
{
	double range = rhi - rlo;

	// XXX not 20, but 1/(1-w) etc etc.
	for (int i=0; i< 25; i++)
	{
		add_to_bin (rlo + range, scale);
		range *= 0.5;
		scale *= weight;
	}

	range = rhi - rlo;
	range *= 0.25;
	for (int i=0; i< 23; i++)
	{
		add_to_bin (rhi - range, scale);
		range *= 0.5;
		scale *= weight;
	}
}

void range (double range, double weight, double scale)
{
	double rlo;
	double delt = 1.0 / range;
	for (rlo = 0.0; rlo < 1.0; rlo += delt)
	{
		hseq (rlo, rlo+delt, weight, scale);
	}
}

void taki(double veight, double weight)
{
	int k;
	double scale = 1.0;
	double tk = 1.0;
	for (k=0; k<16; k++)
	{
		range (tk, weight, scale);
		scale *= veight;
		tk *= 2.0;
	}
}


main(int argc, char *argv[])
{
	int i;

	if (argc < 4)
	{
		fprintf (stderr, "Usage: %s <nbins> <veight> <weight>\n", argv[0]);
		exit (1);
	}
	int nb = atoi (argv[1]);
	double veight = atof (argv[2]);
	double weight = atof (argv[3]);

	alloc_bins (nb);
	printf ("# veight=%g weight = %g\n", veight, weight);
	taki (veight, weight);
	prt_bins();
}

