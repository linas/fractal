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

// The finite-length expander function
double finite_pdr(char* bitseq, int nbits, double Kay, double y)
{
	double acc = y;
	for (int i=0; i<nbits; i++)
	{
		acc *= 1.0 / (2.0*Kay);
		if (bitseq[nbits-i-1])
		{
			acc += 0.5;
		}
	}
	return acc;
}

// The finite-length expander function
// Alternate implementation, for bug-checking
double finite_pdr_alt(char* bitseq, int nbits, double Kay, double y)
{
	double acc = 0.0;
	double tk = 1.0;
	for (int i=0; i<nbits; i++)
	{
		if (bitseq[i]) acc += tk;
		tk *= 1.0 / (2.0*Kay);
	}
	acc = 0.5 * acc + y/tk;

	return acc;
}

// The step function.
int step(char* bitseq, int nbits, double Kay, double y)
{
	if (0 == nbits)
	{
		// if ...
	}
	return 0;
}

// The tree function
double tree_fun(double x, double Kay, double y)
{
	// Decompose x into a bit sequence.
	char bitseq[50];
	for (int i=0; i<50; i++)
	{
		if (0.5 <= x)
		{
			x -= 0.5;
			bitseq[i] = 1;
		}
		else bitseq[i] = 0;
		x *= 2.0;
	}

	double a = finite_pdr(bitseq, 3, Kay, y);
	double b = finite_pdr_alt(bitseq, 3, Kay, y);

	// Construct the tree function from this sequence
	return a-b;
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
