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
mpf_t one;

void do_init(int nbits)
{
	gmp_randinit_default(rstate);
	mpf_set_default_prec(nbits);

	mpf_init(one);
	mpf_init(half);
	mpf_set_ui(one, 1);
	mpf_div_ui(half, one, 2);
}

// Given a double-precision value x, this will create a random
// bit-sequence that is nbits long, with the top 50 bits being
// those taken from the double-precision value x, and the rest
// randomly generated.  This is meant to provide a uniform sampling
// on the unit interval; equivalently, uniform sampling on the
// product space.
void make_random_bitsequence(mpf_t& val, double x, int nbits)
{
	mpf_t tail;
	mpf_init(tail);

	mpf_set_d(val, x);
	mpf_urandomb(tail, rstate, nbits);

	// Keep the top 6 decimal digits of x
	mpf_div_ui(tail, tail, 1000000);
	mpf_add(val, val, tail);
}

void bern(mpf_t& ex, mpf_t TwoKay)
{
	if (0 <= mpf_cmp(ex, half))
	{
		mpf_sub(ex, ex, half);
	}
	mpf_mul(ex, ex, TwoKay);
}

void tent(mpf_t& ex, mpf_t TwoKay)
{
	if (0 <= mpf_cmp(ex, half))
	{
		mpf_ui_sub(ex, 1, ex);
	}
	mpf_mul(ex, ex, TwoKay);
}

void tarp(mpf_t& ex, mpf_t TwoKay)
{
	mpf_mul(ex, ex, TwoKay);
	if (0 <= mpf_cmp(ex, one))
	{
		mpf_sub(ex, TwoKay, ex);
	}
}

void feig(mpf_t& ex, mpf_t four_Kay, mpf_t& scratch)
{
	// K *= 4.0; -- already pased in
	// return 4 * K * x * (1.0 - x);
	mpf_ui_sub(scratch, 1, ex);
	mpf_mul(ex, ex, scratch);
	mpf_mul(ex, ex, four_Kay);
}

void bifurcation_diagram (float *array,
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

	mpf_t Kay, ex, scratch;
	mpf_init(Kay);
	mpf_init(ex);
	mpf_init(scratch);

#define FULLSCALE
#ifdef FULLSCALE
	K = 0.5* (1.0+K);
#endif
#ifdef FEIGEN
	mpf_t four_Kay;
	mpf_init(four_Kay);

	// hack for feigenbaum:
	K = 1.0 - 0.25* (1.0-K);
#endif

	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	int cnt=0;

	for (int j=0; j<itermax; j++)
	{
		double jit = rand();
		jit /= RAND_MAX;
		jit /= array_size;
#ifdef FEIGEN
		jit *= 0.35; /// argh hack for feigenbaum.
#endif

		// The incoming K is always for the top edge of the pixel.
		// The minus sign on the jitter drives it downwards into the pixel.
		make_random_bitsequence(Kay, 2.0*(K-jit), NBITS);
#ifdef FEIGEN
		mpf_mul_ui(four_Kay, Kay, 2);
#endif

		double t = rand();
		t /= RAND_MAX;
		double x = t;
		make_random_bitsequence(ex, x, NBITS);

		/* OK, now start iterating the map */
		for (int iter=0; iter < NBITS; iter++)
		{
			// bern(ex, Kay);
			// tent(ex, Kay);
			tarp(ex, Kay);
			// feig(ex, four_Kay, scratch);

			// Excluding the settle time, together with the support
			// normalization really brings out the details and the color
			// in the tent map. Note that 80 = 2 percent of the 4K
			int settle_time = 80;
			if (settle_time < iter)
			{
				double x = mpf_get_d(ex);
// #define SIDESCALE
#ifdef SIDESCALE
				x /= K;
#endif
// #define SIDETENT
#ifdef SIDETENT
				x = (x - 2.0*K*(1.0-K)) / (K*(2.0*K - 1.0));
#endif
				double en = array_size * (x-floor(x));
				int n = en;
				if (0 > n) n = 0;
				if (n >= array_size) n = array_size-1;
				array[n] += 1.0;
				cnt ++;
			}
		}
	}

#ifdef SUPPORT
	// Measure the support on the row: the total number of pixels
	// that had something in them.
	double support = 0;
	for (int j=0; j<array_size; j++)
		if (0.5 < array[j]) support += 1.0;
	support /= array_size;

	// Bring out the colors
	for (int j=0; j<array_size; j++)
		array[j] *= support;
#endif

// #define SQUARE_INTEGRABLE
#ifdef SQUARE_INTEGRABLE
	// square-integrable norm
	double norm = 0.0;
	for (int j=0; j<array_size; j++)
		norm += array[j] * array[j];

	norm /= array_size;
	norm = 1.0 / sqrt(norm);

	for (int j=0; j<array_size; j++)
		array[j] *= norm;
#endif

#define LP_ONE_NORM
#ifdef LP_ONE_NORM
	// lp_norm for p=1. Interesting but ...
	for (int j=0; j<array_size; j++)
		array[j] *= ((double) array_size) / ((double) cnt);
#endif

#ifdef TENTY
	// Hackery for tent map
	for (int j=0; j<array_size; j++)
		array[j] *= 2.0*(K-0.5);
#endif
#ifdef SIDETENT
	// Hackery for tent map
	for (int j=0; j<array_size; j++)
		// array[j] *= 2.0*(K-0.5);
		array[j] *= 2.0*(2.0-K)*(K-0.5); // meaningless but looks nice
		// array[j] *= K;
#endif

}

DECL_MAKE_BIFUR(bifurcation_diagram)

#if MIDPOINT_GAMES
// Verify the distruibution of the midpoint, so somthing like that ...
static void midpoint_diagram (float *array,
                              int array_size,
                              double x_center,
                              double x_width,
                              double K,
                              int itermax,
                              double omega)
{
	static bool init = false;
	if (not init)
	{
		init = true;
		do_init(itermax);
	}

	mpf_t Kay, ex, twoK;
	mpf_init(Kay);
	mpf_init(ex);
	mpf_init(twoK);

	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	int cnt=0;

	// Compute
	// unsigned long int deno = 128*81*125*49*121*13*17*19*23*29;
	// unsigned long int deno = 64*27*125*49*11*13;
	unsigned long int deno = 32*9*25;
	unsigned long int num = K * ((double) deno);
	while (0 == num%2 and 0 == deno%2) { num /=2; deno /=2; }
	while (0 == num%3 and 0 == deno%3) { num /=3; deno /=3; }
	while (0 == num%5 and 0 == deno%5) { num /=5; deno /=5; }
	printf("#\n# K = %lu / %lu = %g\n#\n", num, deno, K);
	mpf_set_ui(Kay, num);
	mpf_div_ui(Kay, Kay, deno);
	mpf_mul_ui(twoK, Kay, 2);

	mpf_mul(ex, Kay, twoK);
	mpf_sub(ex, ex, Kay);

	/* OK, now start iterating the map */
	for (int iter=0; iter < itermax; iter++)
	{
		double x = mpf_get_d(ex);
		double en = array_size * (x-floor(x));
		int n = en;
		if (0 > n) n = 0;
		if (n >= array_size) n = array_size-1;
		array[n] += 1.0;
		cnt ++;

		if (0 <= mpf_cmp(ex, half))
		{
			mpf_sub(ex, ex, half);
		}
		mpf_mul(ex, ex, twoK);
	}

	// lp_norm for p=1. Interesting but ...
	for (int j=0; j<array_size; j++)
		array[j] *= ((double) array_size) / ((double) cnt);
}

// Compute eigenfunction, recursively.
double reig(double x, double K, int niter)
{
	if (K < x) return 0.0;
	if (niter < 0)
	{
		// Approximate by a constant.
		return 1.0 / K;
	}

	double tkay = 2.0*K;
	double sum = reig(x/tkay, K, niter-1);
	sum += reig(0.5 + x/tkay, K, niter-1);
	sum /= tkay;
	return sum;
}

int main (int argc, char * argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K itermax\n", argv[0]);
		exit(-1);
	}
	double Kay = atof(argv[1]);
	double itermax = atof(argv[2]);

#define ARRSZ 203
	float birr[ARRSZ];
	float arr[ARRSZ];

	bifurcation_diagram (birr, ARRSZ, 0.0, 0.0, Kay, 1200, 0.0);
	fprintf(stderr, "Done with bifur\n");
	midpoint_diagram (arr, ARRSZ, 0.0, 0.0, Kay, itermax, 0.0);
	fprintf(stderr, "Done with mid\n");

	for (int i=0; i<ARRSZ; i++)
	{
		double x = (((double) i) + 0.5) / ((double) ARRSZ);
		double eig = reig(x, Kay, 35);
		printf("%d	%g	%g	%g	%g\n", i, x, arr[i], birr[i], eig);
	}
}
#endif
