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

// #include "brat.h"

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

void feig(mpf_t& ex, mpf_t four_Kay, mpf_t& scratch)
{
	// K *= 4.0; -- already pased in
	// return 4 * K * x * (1.0 - x);
	mpf_ui_sub(scratch, 1, ex);
	mpf_mul(ex, ex, scratch);
	mpf_mul(ex, ex, four_Kay);
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

	mpf_t Kay, ex, scratch, four_Kay;
	mpf_init(Kay);
	mpf_init(ex);
	mpf_init(scratch);
	mpf_init(four_Kay);

#ifdef FEIGEN
	// hack for feigenbaum:
	K = 1.0 - 0.25* (1.0-K);
#endif

	/* clear out the row */
	for (int j=0; j<array_size; j++) array[j] = 0.0;

	int cnt=0;

	for (int j=0; j<itermax; j++)
	{
#ifdef RANDO
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
#endif // RANDO

#define MIDPOINT
#ifdef MIDPOINT
		// Compute 
		// unsigned long int deno = 128*81*125*49*121*13*17*19*23*29;
		// unsigned long int deno = 64*27*125*49*11*13;
		unsigned long int deno = 32*9*25;
		unsigned long int num = K * ((double) deno);
		while (0 == num%2) { num /=2; deno /=2; }
		while (0 == num%3) { num /=3; deno /=3; }
		while (0 == num%5) { num /=5; deno /=5; }
		printf("#\n# K = %lu / %lu = %g\n#\n", num, deno, K);
		mpf_set_ui(Kay, num);
		mpf_div_ui(Kay, Kay, deno);

		mpf_mul(ex, Kay, Kay);
		mpf_mul_ui(ex, ex, 2);
		mpf_sub(ex, ex, Kay);
#endif

		/* OK, now start iterating the map */
		for (int iter=0; iter < NBITS; iter++)
		{
			bern(ex, Kay);
			// tent(ex, Kay);
			// feig(ex, four_Kay, scratch);

#ifdef SETTLE
			// Excluding the settle time, together with the support
			// normalization really brings out the details and the color
			// in the tent map. Note that 80 = 2 percent of the 4K
			int settle_time = 80;
#endif
int settle_time = 0;
			if (settle_time < iter)
			{
				double x = mpf_get_d(ex);
				double en = array_size * (x-floor(x));
				int n = en;
				if (0 > n) n = 0;
				if (n >= array_size) n = array_size-1;
				array[n] += 1.0;
				cnt ++;
			}
		}
	}

	// Measure the support on the row: the total number of pixels
	// that had something in them.
	double support = 0;
	for (int j=0; j<array_size; j++)
		if (0.5 < array[j]) support += 1.0;
	support /= array_size;

	// Bring out the colors
	for (int j=0; j<array_size; j++)
		array[j] *= support;

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
