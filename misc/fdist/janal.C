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
	position -= floor (position);

	// position*nbins then truncate and bincount.
	position *= nbins;
	int i = (int) position;
	bin[i] += amount;
}

void hseq (double rlo, double rhi, double weight, double scale)
{
	double range = rhi - rlo;

	// XXX not 20, but 1/(1-w) etc etc.
	for (int i=0; i< 20; i++)
	{
		add_to_bin (rlo + range, scale);
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

void taki(double weight, double veight)
{
	int k;
	double scale = 1.0;
	double tk = 1.0;
	for (k=0; k<10; k++)
	{
		range (tk, weight, scale);
		scale *= veight;
		tk *= 2.0;
	}
}


main(int argc, char *argv[])
{
	int i;

	if (argc < 2)
	{
		fprintf (stderr, "Usage: %s <nbins>\n", argv[0]);
		exit (1);
	}
	int nb = atoi (argv[1]);

	alloc_bins (nb);
	taki (0.5, 0.5);
	prt_bins();
}

