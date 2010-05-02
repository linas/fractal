/*
 * Eigenfunctions of the dyadic sawtooth, based on the 
 * Walsh functions.
 *
 * Linas Vepstas May 2010
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Square wave function function, returns +1 if x< 1/2 and returns -1 if x> 1/2
 * Repeats with period 1.
 */
void step_1(mpf_t result, mpf_t x)
{
	mpz_t intpart;
	mpf_t ex;
	static int is_init=0;
	static mpf_t half;

	if (0 == is_init)
	{
		is_init = 1;
		mpf_init(half);
		mpf_set_ui(half, 1);
		mpf_div_ui(half, half, 2);
	}

	mpf_init(ex);
	mpf_set(ex, x);
	mpz_init(intpart);
	mpz_set_f(intpart, ex);
	mpf_set_z(result, intpart);
	mpf_sub(result, ex, result);

	if (0 < mpf_cmp(half, result))
	{
		mpf_set_ui(result, 1);
		return;
	}
	
	mpf_set_si(result, -1);
}

/**
 * returns square wave of frequency 2^(n-1)
 * n must be less than 32/64
 * x must be in range of zero, one
 */
void step_n(mpf_t result, mpf_t x, int n)
{
	unsigned long tp;

	if (1 == n)
	{
		step_1(result, x);
		return;
	}

	tp = 1<<(n-1);
	mpf_set_ui(result, tp);
	mpf_mul(result, result, x); 
	step_1(result, result);
}

/**
 * Implement the n'th walsh function
 * n must be less that 2^32 or 2^64
 */
void walsh(mpf_t result, mpf_t x, unsigned long n)
{
	int i;
	mpf_t term;
	mpf_init(term);

	mpf_set_ui(result, 1);
	for (i=1; i<=32; i++)
	{
		if (n & 0x1)
		{
			step_n(term, x, i);
			mpf_mul(result, result, term);
		}

		n >>= 1;
	}
}

int main (int argc, char * argv[])
{
	mpf_t x, y, step;
	double x_f, y_f;
	int i, npts;
	int prec, nbits;

	if (2 > argc)
	{
		fprintf(stderr, "Usage: %s <decimal-precision>\n", argv[0]);
		exit(1);
	}

	/* prec is decimal-places of precision */
	prec = 50;
	prec = atoi(argv[1]);

	/* Set the precision (number of binary bits) */
	nbits = 3.3*prec;
	mpf_set_default_prec (nbits);

	mpf_init(x);
	mpf_init(y);
	mpf_init(step);

	npts = 100;
	mpf_set_ui(step, 1);
	mpf_div_ui(step, step, npts);

	mpf_set_ui(x, 0);
	for (i=0; i<npts; i++)
	{
		// step_1(y, x);
		// step_n(y, x, 3);
		walsh(y, x, 6);

		x_f = mpf_get_d(x);
		y_f = mpf_get_d(y);

		printf("%d	%f	%g\n", i, x_f, y_f);

		mpf_add(x, x, step);
	}

	return 0;
}
