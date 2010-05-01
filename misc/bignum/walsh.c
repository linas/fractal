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
	static int is_init=0, mpf_t half;
	if (0 == is_init)
	{
		is_init = 1;
		mpf_init(half);
		mpf_set_ui (half, 1);
		mpf_div_ui (half, half, 2);
	}
}

main ()
{
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

}
