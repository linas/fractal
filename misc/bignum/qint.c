
/*
 * integral of question mark measure.
 *
 * Linas March 2010
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mp-quest.h"

int main (int argc, char * argv[])
{
	int prec, nbits;
	mpf_t x, q;
	double y;

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
	mpf_init(q);
	mpf_set_ui (x, 1);
	mpf_div_ui (x, x, 3);

	question_mark(q, x, prec);

	y = mpf_get_d(q);

	printf("its %g\n", y);
	
	return 0;
}
