
/*
 * ftag.c
 *
 * Rebuild the Farey dist as a Cantor set
 * This time, as a log of Takagi thing.
 * Use a stupa energy potential, compare and contrast to
 * the Kac potential
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

/* A stupa is a step pyramid */
double stupa (double x)
{
	int n;

	x -= floor(x);
	if (0.5 < x) x = 1.0-x;
	if (0.0 == x) return 1.0e16;
	if (0.5 == x) return 1.0e16;
	
	n = 0;
	while(1)
	{
		if (0.25 < x) break;
		n++;
		x *= 2.0;
	}
	return (double) n;
}

/* The Kac potential */
double kac (double x)
{
	x -= floor(x);
	if (x < 0.5) return 2.0*x-0.5;
	return 1.5 - 2.0*x;
}

double taga(double x, double veight)
{
	int k;

	double acc = 0.0;
	double tk = 1.0;
	double tw = 1.0;
	for (k=0; k<30; k++)
	{
		// double f = stupa(tk*x);
		double f = kac(tk*x);
		acc += tw * f;
		tw *= veight;
		tk *= 2.0;
	}

	return acc;
}

double product(double x, double veight, double weight)
{
	double lg = taga(x, veight);
	lg *= log (weight);
	lg = exp(lg);
	return lg;
}

void prt_graph(int npts, double veight, double weight)
{
	int i;

	double *bin = (double *) malloc(npts*sizeof(double));
	for (i=0; i<npts; i++)
	{
		bin[i] = 0.0;
	}

	ContinuedFraction f;

	/* Compute the integral of the distribution, in "gral" */
	/* run loop once, to get the total normalization */
	double gral = 0.0;
	double delta = 1.0 / ((double) npts);
	for (i=0; i<npts; i++)
	{
		double x = ((double) i) / ((double) npts);

#if NO_JITTER
   	f.SetRatio (i, npts);
   	double far = f.ToFarey (); 
		double p = product (far, veight, weight);
		gral += p * delta;
		bin[i] += p;
#endif
#define JITTER
#ifdef JITTER
		for (int j=0; j<1120; j++)
		{
			double off = delta * rand() / ((double) RAND_MAX);
   		f.SetReal (x+off);
   		double far = f.ToFarey (); 
			double p = product (far, veight, weight);
			gral += p * delta;
			bin[i] += p;
		}
#endif
	}

	delta /= gral;
	gral = 0.0;
	double entropy = 0.0;
	for (i=0; i<npts; i++)
	{
		double x = ((double) i) / ((double) npts);
		double p = bin[i];
		gral += p * delta;
		if (p != 0.0) entropy -= p * delta * log(p);

   	f.SetRatio (i, npts);
   	double far = f.ToFarey (); 

		printf ("%d	%8.6g	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, p, gral, far, entropy);
	}
	printf ("# Total entropy = %g (log2=%g)\n", entropy, log(2.0)*entropy);
}

int main (int argc, char * argv[])
{
	if (argc < 4)
	{
		fprintf (stderr, "Usage: %s <nbins> <veight> <weight>\n", argv[0]);
		exit (1);
	}
	int nbins = atoi (argv[1]);
	double veight = atof (argv[2]);
	double weight = atof (argv[3]);

	prt_graph (nbins, veight, weight);

	return 0;
}
