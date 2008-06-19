/*
 * Compute matrix elements of fourier of minkowski question mark 
 * by means of brute-force integration. This seems to be a viable
 * quick-n-dirty way of getting these.
 *
 * Linas June 2008
 */

#include <math.h>

static int arrsz;
static long double *sin_arr;
static long double *cos_arr;
static long double delta;

void init_sine_cosine_arrays(int size)
{
	arrsz = size;
	sin_arr = (long double *) malloc(arrsz * sizeof(long double));
	cos_arr = (long double *) malloc(arrsz * sizeof(long double));

	delta = 2.0L*M_PIl;
	delta /= (long double) arrsz;
	for (int i=0; i<arrsz; i++)
	{
	}

}
