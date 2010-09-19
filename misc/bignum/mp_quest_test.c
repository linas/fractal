
/*
 * Small unit-test module for the MP minkowski-question-mark function
 *
 * Copyright (C) 2010 Linas Vepstas <linasvepstas@gmail.com>
 */

#include <gmp.h>
#include <stdio.h>
#include "mp-quest.h"

int
main()
{
	int prec, nbits;
	mpf_t x, y;
	double x_f, y_f;

	prec = 90;
	/* Set the precision (number of binary bits) */
	nbits = 3.3*prec;
	mpf_set_default_prec (nbits);

	mpf_init(x);
	mpf_init(y);

	mpf_set_ui(x, 959);
	mpf_div_ui(x, x, 1233);

	question_mark(y, x, prec);

	x_f = mpf_get_d(x);
	y_f = mpf_get_d(y);

	printf("its %f %f\n", x_f, y_f);

	return 0;
}
