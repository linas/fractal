
/*
 * Conveience to graph the a_k
 *
 * Linas Vepstas October 2008
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double a_k (double x, int k)
{
	int i;
	for (i=0; i<k; i++)
	{
		if (x < 0.5) x = x/(1.0-x);
		else x = (2.0*x-1)/x; 
		// else x = (1.0-x)/x; 
	}
	return x;
}

void graph(int npts)
{
	int i;

	for (i=0; i<npts; i++)
	{
		double x = (double) i / ((double) npts);

		printf ("%g	%g %g	%g	%g\n", x, a_k(x,0), a_k(x,1), a_k(x,2), a_k(x,3));
		fflush (stdout);
	}
}


main(int argc, char *argv[])
{
	int i;

	if (argc <2)
	{
		fprintf (stderr, "Usage: %s <nbins> \n", argv[0]);
		exit (1);
	}
	int nbins = atoi (argv[1]);

	graph (nbins);
}
