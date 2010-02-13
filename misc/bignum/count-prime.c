/*
 * count-prime.c
 *
 * A prime-counting conjecture.
 *
 * Copyright (C) 2010 Linas Vepstas
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 *
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mp-quest.h"

int main (int argc, char * argv[])
{
	unsigned long int n;
	mpf_t qi, x, acc;
	int i, prec, nbits, step, twostep, prt;
	int mmax;

	mpf_t *qinv;

	if (3 > argc)
	{
		fprintf(stderr, "Usage: %s <decimal-precision> <1/step>\n", argv[0]);
		exit(1);
	}

	/* prec is decimal-places of precision */
	/* step is the integral step size */
	prec = 50;
	prec = atoi(argv[1]);
	step = atoi(argv[2]);

	/* mmax is how far we weill go */
	/* 2^10 = 10^3 approx fo 2^80 = 10^24 */
	mmax = 80;

	printf("#\n# Prime-counting-conjecture via question mark\n#\n");
	printf("#\n# Decimal precision = %d\n#\n", prec);
	printf("#\n# step size = 1 / %d\n#\n", step);

	/* Set the precision (number of binary bits) */
	nbits = 3.3*prec;
	mpf_set_default_prec (nbits);

	/* Alloc array that will hold questioin-mark values */
	qinv = (mpf_t *) malloc (step * sizeof(mpf_t));

	mpf_init(x);
	mpf_init(qi);

	/* compute 1 / ( ?^{-1}(1/x) ) for 1 =< x < 2 */
	for(i=0; i<step; i++)
	{
		mpf_init(qinv[i]);

		mpf_set_ui (x, step);
		mpf_div_ui (x, x, step+i);
		question_inverse(qi, x, prec);
		mpf_ui_div (qinv[i], 1, qi);
	}


	mpf_init(acc);

	mpf_set_ui(acc, 0);

	twostep = 2*step;
	n = twostep;
	prt = step / 10;
	while(1)
	{
		mpf_set_ui (x, twostep);
		mpf_div_ui (x, x, n);

		question_inverse(qi, x, prec);
		mpf_add (acc, acc, qi);

		if (n%prt == 0)
		{
			double xd = mpf_get_d(x);
			xd = 2.0 / xd;
			double facc = mpf_get_d(acc);
			facc /= (double) step;
			facc /= log(2.0);
			printf("%lu	%f	%f\n", n, xd, facc);
			fflush(stdout);

			if (n > 10*prt) prt *= 10;
		}
		n++;
	}
	
	return 0;
}

/* =============================== END OF FILE =========================== */
