/*
 * rand-measure.C
 *
 * Distribution of the random numbers
 * aka multi-fractal measures
 *
 * Linas January 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "gcf.h"

double randoid(int n)
{
#define NVAL 2265186
	static double array[NVAL];
	static int inited=0;
	if (!inited)
	{
		int i;
		inited = 1;
		srand (99);
		for (i=0; i<NVAL; i++)
		{
			// array[i] = rand() & 0x1;
			array[i] = ((double) rand()) / ((double)RAND_MAX);
		}
	}

	if (NVAL<= n) return 0;
	double val = array[n];

#if 0
	val *= n;
	val *= n;
	val *= n;
#endif

	return val;
}

main(int argc, char *argv[])
{
	int i;

	if (argc <2)
	{
		fprintf (stderr, "Usage: %s <nbins>\n", argv[0]);
		exit (1);
	}
	int nbins = atoi (argv[1]);

	double z = 0.0;
	for (i=0; i<nbins; i++)
	{
		double x = ((double) i) / ((double) nbins);

		int g = gcf32 (i, nbins);
		int p = i/g;
		int q = nbins/g;

		double y = randoid (p*q);
		z += y / ((double) nbins);

		printf ("%d	%g	%g	%g\n", i,x ,y,z);
	}
}

