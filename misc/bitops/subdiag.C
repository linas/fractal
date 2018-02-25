/*
 * subdiag.C
 *
 * Goofy color plot of the subdiagonal matrix elements <n+1|L|n>  
 *
 * Is not very interesting -- some structure for small n, and mixed,
 * random mess for the higher n.
 *
 * February 2018
 */

#include <math.h>
#include <stdio.h>
#include <pthread.h>

#include "brat.h"

#define NOMAIN
#include "psi.c"
#include "psibig.c"

static void subdiagonal(float *array,
                        int array_size,
                        double x_center,
                        double x_width,
                        double K,
                        int itermax,
                        double omega)
{
	static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	pthread_mutex_lock(&mutex);

	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	double Kay = 0.5 + 0.5*K;

	// find_midpoints(Kay);
	big_midpoints(Kay, 400, midpoints, MAXN);
	sequence_midpoints(Kay, MAXN);

	for (int j=0; j<array_size; j++)
		array[j] = hess(Kay, j+1, j);

	// printf("duuude  K=%g ar=%g\n", Kay, hess(Kay, 3, 2));

	pthread_mutex_unlock(&mutex);
}

DECL_MAKE_BIFUR(subdiagonal)
