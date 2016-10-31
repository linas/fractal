/*
 * Circular cut/integral.
 *
 * Oct 2016
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <mp-arith.h>
#include <mp-trig.h>
#include <moebius.h>

#include "genfunc.h"

int main()
{


	double rad = 6000.0;
	int nsteps = 2000;

	printf("#\n# radius=%g\n#\n", rad);

	int nprec = 85;

	mpf_set_default_prec(nprec * 3.321+ 50);

	cpx_t val, z; cpx_init(val); cpx_init(z);
	mpf_t mag; mpf_init(mag);
	mpf_t delta; mpf_init(delta);
	mpf_t theta; mpf_init(theta);
	mpf_t radius; mpf_init(radius);
	mpf_t zero; mpf_init(zero);

	mpf_t sum; mpf_init(sum);

	double delt = 2.0*M_PI / ((double)nsteps);
	double thet = 0.0;

	mpf_set_ui(zero, 0);
	mpf_set_ui(theta, 0);
	mpf_set_ui(sum, 0);
	mpf_set_d(delta, delt);
	mpf_set_d(radius, rad);

	for (int i=0; ; i++)
	{
		// Compute z = r exp(itheta)
		cpx_set_mpf(z, zero, theta);
		cpx_exp(z, z, nprec);
		cpx_times_mpf(z, z, radius);

		// cpx_exponential_genfunc_mpf(val, z, nprec, parti_z_mpf);
		// cpx_exponential_genfunc(val, z, nprec, big_omega);
		cpx_exponential_genfunc(val, z, nprec, divisor);

		// Sum the magnitudes
		cpx_abs(mag, val);
		mpf_add(sum, sum, mag);

		double th = mpf_get_d(theta);
		double mg = mpf_get_d(mag);
		double sm = mpf_get_d(sum);
		sm *= delt;

		printf("%d	%g	%g	%g\n", i, th, mg, sm);

		mpf_add(theta, theta, delta);
		thet += delt;
		if (2.0*M_PI < thet) break;
	}
}
