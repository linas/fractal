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
	double qid;
	mpf_t qi, x, acc;
	int prec, nbits, step, twostep, prt;

	if (3 > argc)
	{
	`	fprintf(stderr, "Usage: %s <decimal-precision> <1/step>\n", argv[0]);
		exit(1);
	}

	prec = 50;
	prec = atoi(argv[1]);
	step = atoi(argv[1]);

	printf("#\n# Prime-counting-conjecture via question mark\n#\n");
	printf("#\n# Decimal precision = %d\n", prec);

	/* Set the precision (number of binary bits) */
	nbits = 3.3*prec;
	mpf_set_default_prec (nbits);

	mpf_init(x);
	mpf_init(qi);
	mpf_init(acc);

	mpf_set_ui(acc, 0);

	n = step;
	twostep = 2*step;
	prt = step / 10;
	while(1)
	{
		mpf_set_ui (x, twostep);
		mpf_div_ui (x, x, n);

		question_inverse(qi, x, prec);
		mpf_add (acc, qi, qi);

		if (n%prt == 0)
		{
			double xd = mpf_get_d(x);
			double facc = mpf_get_d(acc);
			printf("%ud	%f	%f\n", n, xd, facc);

			if (n > 10*prt) prt *= 10;
		}
		n++;
	}
	
	return 0;
}

/* =============================== END OF FILE =========================== */
