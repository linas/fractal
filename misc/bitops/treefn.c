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
	acc = 0.5 * acc + y * tk;

	return acc;
}

// The step function.
int step(char* bitseq, int nbits, double Kay, double y)
{
	double ve = finite_pdr(bitseq, nbits, Kay, y);
	if (ve <= Kay) return 1;
	return 0;
}

// Compute the bit sequence of x.
void to_bit_sequence(char* bitseq, double x)
{
	// Decompose x into a bit sequence.
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
}

// The gamma function
double gamma_fun(double x, int nbits, double Kay, double y)
{
	// Decompose x into a bit sequence.
	char bitseq[50];
	to_bit_sequence(bitseq, x);

	return finite_pdr(bitseq, nbits, Kay, y);
}

// The tree function
double tree_fun(double x, double Kay, double y)
{
	// Decompose x into a bit sequence.
	char bitseq[50];
	to_bit_sequence(bitseq, x);

	// Construct the tree function from this sequence
	for (int i=0; i<50; i++)
	{
		if (0 == step(bitseq, i, Kay, y)) return 0.0;
	}

	return 1.0;
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

	int npts = 1603;
	for (int i=0; i<npts; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) npts);
		printf("%d	%g", i, x);
#if SLICES
		for (why=0.0; why<=1.0; why += 0.1)
		{
			double t = tree_fun(x, Kay, why);
			printf("	%g", t);
		}
		printf("\n");
#endif
		double t = tree_fun(x, Kay, why);
		printf("	%g", t);
		for (int n=1; n<11; n++)
		{
			double g = gamma_fun(x, n, Kay, why);
			printf("	%g", g);
		}
		printf("\n");
	}
}
