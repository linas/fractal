/*
 * bigbern.C
 *
 * Bignum altered Bernoulli map
 * Dec 2017
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/*-------------------------------------------------------------------*/
/*
 */

gmp_randstate_t rstate;
mpf_t half;

void do_init(int nbits)
{
	gmp_randinit_default(rstate);
	mpf_set_default_prec(nbits);

	mpf_init(half);
	mpf_set_ui(half, 1);
	mpf_div_ui(half, half, 2);
}

void make_mpf(mpf_t& val, double x, int nbits)
{
	mpf_t tail;
	mpf_init(tail);

	mpf_set_d(val, x);
	mpf_urandomb(tail, rstate, nbits);

	// Keep the top 6 decimal digits of x
	mpf_div_ui(tail, tail, 1000000);
	mpf_add(val, val, tail);
}

void bern(mpf_t& ex, mpf_t Kay)
{
	if (0 <= mpf_cmp(ex, half))
	{
		mpf_sub(ex, ex, half);
	}
	mpf_mul(ex, ex, Kay);
}

void tent(mpf_t& ex, mpf_t Kay)
{
	if (0 <= mpf_cmp(ex, half))
	{
		mpf_ui_sub(ex, 1, ex);
	}
	mpf_mul(ex, ex, Kay);
}

double feig(double x, double K)
{
	K *= 4.0;
	return K * x * (1.0 - x);
}

static void bifurcation_diagram (float *array,
                                 int array_size,
                                 double x_center,
                                 double x_width,
                                 double K,
                                 int itermax,
                                 double omega)
{
#define NBITS 4000

	static bool init = false;
	if (not init)
	{
		init = true;
		do_init(NBITS);
	}

	mpf_t Kay, ex;
	mpf_init(Kay);
	mpf_init(ex);

	make_mpf(Kay, 2.0*K, NBITS);

	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	int cnt=0;

	for (int j=0; j<itermax; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;
		make_mpf(ex, x, NBITS);

		/* OK, now start iterating the benoulli map */
		for (int iter=0; iter < NBITS; iter++)
		{
			bern(ex, Kay);
			// tent(ex, Kay);
			// x = feig(x, K);

			x = mpf_get_d(ex);
			double en = array_size * (x-floor(x));
			int n = en;
			if (0 > n) n = 0;
			if (n >= array_size) n = array_size-1;
			array[n] += 1.0;
			cnt ++;
		}
	}

	for (int j=0; j<array_size; j++)
		array[j] *= ((double) array_size) / ((double) cnt);
}

// DECL_MAKE_BIFUR(bifurcation_diagram)

#if 1
int main (int argc, char * argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K itermax\n", argv[0]);
		exit(-1);
	}
	double Kay = atof(argv[1]);
	double itermax = atof(argv[2]);

#define ARRSZ 803
	float arr[ARRSZ];

	bifurcation_diagram (arr, ARRSZ, 0.0, 0.0, Kay, itermax, 0.0);

	double sum = 0;
	for (int i=0; i<ARRSZ; i++)
	{
		double x = (((double) i) + 0.5) / ((double) ARRSZ);
		double rho = arr[i];
		sum += rho / ((double) ARRSZ);
		printf("%d	%g	%g	%g\n", i, x, rho, sum);
	}
}
#endif
