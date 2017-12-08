/*
 * bigbern.C
 *
 * Bignum altered Bernoulli map
 * Dec 2017
 */

#iclude <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/*-------------------------------------------------------------------*/
/*
 */

gmp_randstate_t rstate;

void init()
{
	gmp_randinit_default(rstate);
}

void make_mpf(mpf_t& val, double x, int nbits)
{
	mpz_t digs, tn;
	mpz_init2(tn, nbits);
	mpz_init2(digs, nbits);
	mpz_urandomb(digs, rstate, nbits);
	mpz_ui_pow_ui(tn, 2, nbits);
	
}

double bern(double x, double K)
{
	K *= 2.0;
	if (0.5 <= x)
	{
		return K * (x - 0.5);
	}
	return K*x;
}

double tent(double x, double K)
{
	K *= 2.0;
	if (0.5 <= x)
	{
		return K * (1.0 - x);
	}
	return K*x;
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
	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	int cnt=0;

	for (int j=0; j<itermax; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		/* OK, now start iterating the benoulli map */
		for (int iter=0; iter < 550; iter++)
		{
			// x = bern(x, K);
			// x = noadd(x, K);
			// x = tent(x, K);
			// x = notent(x, K);
			// x = feig(x, K);
			x = nofeig(x, K);

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

DECL_MAKE_BIFUR(bifurcation_diagram)

#if 0
int main ( int argc, char * argv[])
{
	double om = atof(argv[1]);
	double kb = atof(argv[2]);
}
#endif
