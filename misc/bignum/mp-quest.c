/*
 * mp-quest.c
 *
 * High-precison Minkowski Question mark, Stern-Brocot Tree, etc.
 * using the Gnu Multiple-precision library.
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
#include "mp-quest.h"

void question_inverse (mpf_t qinv, const mpf_t x, unsigned int prec)
{
	int i, n;
	unsigned long idx;
	mpf_t mantissa;
	mpz_t bits;

	mpf_init(mantissa);
	mpz_init(bits);

	/* Get the number of binary bits from prec = log_2 10 * prec */
	int nbits = (int) floor (3.322 * prec);
	
   int *bitcnt = (int *) malloc ((nbits+1) * sizeof(int));
	memset (bitcnt, 0, (nbits+1) * sizeof(int));

	mpf_mul_2exp(mantissa, x, nbits);
	mpz_set_f(bits, mantissa);

	/* Count the number of contiguous bits */
	idx = 0;
	n = 0;
	while (n < nbits)
	{
		if (n%2 == 0)
		{
   		idx = mpz_scan0(bits, idx);
		}
		else
		{
   		idx = mpz_scan1(bits, idx);
		}
		if (ULONG_MAX == idx) break;
		bitcnt[n] = idx;
		idx ++;
		n++;
	}


	mpf_set_ui(qinv, 0);

	/* Compute the corresponding continued fraction */
	for (i=0; i<n; i++)
	{
		if (0 == bitcnt[i]) continue;
		mpf_add_ui(qinv, qinv, bitcnt[i]);
		mpf_ui_div(qinv, 1, qinv);
	}
	
   free (bitcnt);
}

/* =============================== END OF FILE =========================== */
