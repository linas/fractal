/*
 * Generating functions
 *
 * Oct 2016
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <mp-arith.h>

#include "genfunc.h"

void parti_z_mpf(mpf_t res, long n)
{
   mpz_t part; mpz_init(part);
   partition_z(part, n);
   mpf_set_z(res, part);
#if 0
static int last=0;
if (last < n) {
printf("duuude parti n=%d bits=%lu\n", n, mpz_sizeinbase(part, 2));
last = n;
}
#endif

   mpz_clear(part);
}

int main()
{
	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);

	int nprec = 285;

	double re_q = 0.0;
	double im_q = 0.0;
	for (re_q = 50; re_q < 700; re_q += 30)
	{
		cpx_set_d(z, re_q, im_q);

		cpx_exponential_genfunc_mpf(sum, z, nprec, parti_z_mpf);
		cpx_abs(val, sum);
		double rv = mpf_get_d(val);

		// printf("r=%g rv =%g\n", re_q, rv);
	}
}
