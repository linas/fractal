
/*
 * ftag.c
 *
 * Rebuild the Farey dist as a Cantor set
 * This time, as a log of Takagi thing.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* A stupa is a step pyramid */
double stupa (double x)
{
	int n;

	x -= floor(x);
	if (0.5 < x) x = 1.0-x;
	if (0.0 == x) return 1.0e16;
	
	n = 0;
	while(1)
	{
		if (0.25 < x) break;
		n++;
		x *= 2.0;
	}
	return (double) n;
}

double taga(double x, double veight)
{
	int k;

	double acc = 0.0;
	double tk = 1.0;
	double tw = 1.0;
	for (k=0; k<20; k++)
	{
		double f = stupa(tk*x);
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

	double gral = 0.0;
	double delta = 1.0 / ((double) npts);
	for (i=0; i<npts; i++)
	{
		double x = ((double) i) / ((double) npts);
		double p = product (x, veight, weight);
		gral += p * delta;

		printf ("%d	%8.6g	%8.6g	%8.6g\n", i, x, p, gral);
	}
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
}
