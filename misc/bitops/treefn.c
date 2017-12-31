/*
 * treefn.c
 *
 * Understand the tree function.
 *
 * Linas Vepstas Dec 2017
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// The tree function
double tree_fun(double x, double Kay, double y)
{
	// Decompose x into a bit sequence.
	char nbits[50];
	for (int i=0; i<50; i++)
	{
		if (0.5 <= x)
		{
			x -= 0.5;
			nbits[i] = 1;
		}
		else nbits[i] = 0;
		x *= 2.0;
	}

	// Construct the tree function from this sequence
	double acc = 0.1;
	for (int i=0; i<50; i++)
	{
		acc *= 1.0 / (2.0*Kay);
		if (nbits[50-i-1])
		{
			acc += 0.5;
		}
	}

	return acc;
}

int main (int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K y\n", argv[0]);
		exit (1);
	}
	double Kay = atof(argv[1]);
	double why = atof(argv[2]);

	int npts = 803;
	for (int i=0; i<=npts; i++)
	{
		double x = ((double) i)/ ((double) npts);
		double t = tree_fun(x, Kay, why);
		printf("%d	%g	%g\n", i, x, t);
	}
}
