/*
 * sidorov-density.C
 * Show the extened density measure in a bifurcation-style graph
 *
 * Linas Vepstas Sept 2020
 */

#include "brat.h"

#include "extended-density.C"

static void ext_measure (float *array,
                         int array_size,
                         double x_center,
                         double x_width,
                         double Kay,
                         int itermax,
                         double param)
{
	int nbits = itermax;

	fprintf(stderr, "Working K=%g nbits=%d\n",  Kay, nbits);

/*
 * Use this for quicker results.
 */
#define MAXDEPTH 6
#define NSAMP 32

/*
#define MAXDEPTH 7
#define NSAMP 64
*/

	double arr[array_size];
	extended_measure(2.0*Kay, MAXDEPTH, NSAMP, arr, array_size, nbits);

	for (int j=0; j<array_size; j++) array[j] = arr[j];
}

DECL_MAKE_BIFUR(ext_measure)

// ================================================================
