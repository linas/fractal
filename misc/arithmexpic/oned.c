/*
 * Step in one direction and print the result.
 *
 * Oct 2016
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <mp-arith.h>
#include <moebius.h>

#include "genfunc.h"

void parti_z_mpf(mpf_t res, long n)
{
   mpz_t part; mpz_init(part);
   partition_z(part, n);
   mpf_set_z(res, part);
#if 0
static int last=0;
if (last < n) {
int bits = mpz_sizeinbase(part, 2);
double ln2t = log(2.0) / log(10.0);
int decs = bits * ln2t;
printf("partiion n=%ld bits=%d prec=%d\n", n, bits, decs);
last = n;
}
#endif

   mpz_clear(part);
}


int main()
{
	int nprec = 85;

	mpf_set_default_prec(nprec * 3.321+ 50);

	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);


	double re_q = 0.0;
	double im_q = 0.0;
	for (re_q = 4500; re_q < 16700; re_q += 100)
	{
		im_q = re_q; // / 3.0;
		cpx_set_d(z, re_q, im_q);

		// cpx_exponential_genfunc_mpf(sum, z, nprec, parti_z_mpf);
		cpx_exponential_genfunc(sum, z, nprec, big_omega);
		cpx_abs(val, sum);
		double rv = mpf_get_d(val);

		// printf("r=%g rv =%g\n", re_q, rv);
	}
}
