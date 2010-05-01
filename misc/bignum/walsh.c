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

void step(mpf_t result, mpf_t x)
{
	static int is_init=0;
	static mpf_t half;

	if (0 == is_init)
	{
		is_init = 1;
		mpf_init(half);
		mpf_set_ui (half, 1);
		mpf_div_ui (half, half, 2);
	}

	if (0 < mpf_cmp(x, half))
	{
		mpf_set_ui(result, 1);
		return;
	}
	
	mpf_set_si(result, -1);
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
		step(y, x);

		x_f = mpf_get_d(x);
		y_f = mpf_get_d(y);

		printf("%d	%f	%g\n", i, x_f, y_f);

		mpf_add(x, x, step);
	}
}
