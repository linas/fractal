/*
 * psidelt.C
 * Verify some form of completeness.
 *
 * Ferbruary 2018
 *
 */
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "brat.h"

#define MAXN 4000
#define NOMAIN
#include "psi.c"
#include "psibig.c"

static void diagonal_diagram (float *array,
                             int array_size,
                             double x_center,
                             double x_width,
                             double row,
                             int itermax,
                             double K)
{
	static bool init=false;
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	if (not init)
	{
		pthread_mutex_lock(&mutex);
		if (not init)
		{
			// find_midpoints(K, MAXN);
			big_midpoints(K, 400, midpoints, MAXN);
			sequence_midpoints(K, MAXN);
			printf("working K=%g\n", K);
			init = true;
		}
		pthread_mutex_unlock(&mutex);
	}

	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	double y = row;
	for (int j=0; j<array_size; j++)
	{
		double x = ((double) j + 0.5) / ((double) array_size);
		double sum = 0.0;
		for (int m=0; m<itermax; m++)
		{
			for (int n=0; n<itermax; n++)
			{
				sum += psi_n(y, K, m) * psi_n(x, K, n);
			}
		}
 		array[j] = sum;
	}
}

DECL_MAKE_BIFUR(diagonal_diagram)
