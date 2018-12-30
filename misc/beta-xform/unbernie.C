/*
 * unbernie.C
 *
 * Altered simplified algorithmic Bernoulli map
 * Dec 2017
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/*-------------------------------------------------------------------*/
/*
 */

double unbmap(double y, double K)
{
	// Iterate on y using mashed Bernoulli, and extract symbol dynamics
	char nbits[50];
	for (int i=0; i<50; i++)
	{
		if (0.5 <= y)
		{
			y -= 0.5;
			nbits[i] = 1;
		}
		else nbits[i] = 0;
		y *= K;
	}

	// Reconstruct x in a mashed bernoulli sequence.
	double acc = 0.1;
	for (int i=0; i<50; i++)
	{
		if (nbits[50-i-1])
		{
			acc += 0.5;
		}
		acc *= 0.5;
	}
	return acc;
}

static void mapping_diagram (float *array,
                             int array_size,
                             double x_center,
                             double x_width,
                             double K,
                             int itermax,
                             double omega)
{
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	int cnt=0;

	for (int j=0; j<itermax; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		x = unbmap(x, 2.0*K);

		x *= 2.0;

		double en = array_size * (x-floor(x));
		int n = en;
		if (0 > n) n = 0;
		if (n >= array_size) n = array_size-1;
		array[n] += 1.0;
		cnt ++;
	}

	for (int j=0; j<array_size; j++)
		array[j] *= ((double) array_size) / ((double) cnt);
}

DECL_MAKE_BIFUR(mapping_diagram)
