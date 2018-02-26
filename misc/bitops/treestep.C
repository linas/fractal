/*
 * treestep.C
 *
 * Tree function visualization
 * Dec 2017
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/*-------------------------------------------------------------------*/

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

/*-------------------------------------------------------------------*/
/*
 */

// #define STEPS
#ifdef STEPS
static void step_diagram (float *array,
                                 int array_size,
                                 double x_center,
                                 double x_width,
                                 double y,
                                 int itermax,
                                 double Kay)
{
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	for (int j=0; j<array_size; j++)
	{
		double x = (((double) j) + 0.5)/ ((double) array_size);
		x -= 0.5;
		x *= x_width;
		x += x_center;
		array[j] = tree_fun(x, Kay, y);
	}
}

DECL_MAKE_BIFUR(step_diagram)
#endif

#ifdef TREE
static void tree_diagram (float *array,
                                 int array_size,
                                 double x_center,
                                 double x_width,
                                 double Kay,
                                 int itermax,
                                 double why)
{
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	for (int j=0; j<array_size; j++)
	{
		double x = (((double) j) + 0.5)/ ((double) array_size);
		x -= 0.5;
		x *= x_width;
		x += x_center;
		array[j] = tree_fun(x, Kay, why);

		for (int i=0; i<2499; i++)
		{
			double t = rand();
			t /= RAND_MAX;
			t -= 0.5;
			t /= (double) array_size;
			t *= x_width;
			double ex = x+t;
			array[j] += tree_fun(ex, Kay, why);
		}
		array[j] /= 2500.0;
	}
}

DECL_MAKE_BIFUR(tree_diagram)
#endif

#define STEP_COLOR
#ifdef STEP_COLOR
static void step_color_diagram (float *array,
                                 int array_size,
                                 double x_center,
                                 double x_width,
                                 double K,
                                 int itermax,
                                 double param)
{
	double Kay = 0.5 + 0.5*K;

	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	for (int j=0; j<array_size; j++)
	{
		double x = (((double) j) + 0.5)/ ((double) array_size);
		x -= 0.5;
		x *= x_width;
		x += x_center;

#define NCOLORS 225
		double y = 0.0;
		for (int n=0; n<NCOLORS; n++)
		{
			y = (((double) n) + 0.5)/ ((double) NCOLORS);
			y = 1.0 - y;
			double rv = tree_fun(x, Kay, y);
			if (0.5 < rv) break;
		}
		array[j] = y;
	}
}

DECL_MAKE_BIFUR(step_color_diagram)
#endif
