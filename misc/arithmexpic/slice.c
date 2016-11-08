/*
 * Circular cut/integral.
 * Consider a circle centered on the origin. Graph values
 * of the function on that circle.
 *
 * Oct 2016
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mp-arith.h>
#include <mp-trig.h>
#include <moebius.h>
#include <totient.h>

#include "genfunc.h"

int main(int argc, char * argv[])
{
	if (4 != argc)
	{
		fprintf(stderr, "Usage: %s <name> <nsteps> <radius>\n", argv[0]);
		exit (1);
	}

	char * name = argv[1];
	int nsteps = atoi(argv[2]);
	double rad = atof(argv[3]);

	printf("#\n# radius=%g\n", rad);
	printf("# nsteps=%d\n#\n", nsteps);

	printf("#\n# name = %s\n#\n", name);

	printf("# Legend:\n");
	printf("# i, theta, magnitude, sum-of-magnitude, sum-of-phase-wrap, phase\n");
	printf("#\n");

	long (*func)(long) = NULL;
	if (0 == strcmp(name, "totient")) func = totient_phi;
	if (0 == strcmp(name, "divisor")) func = divisor;
	if (0 == strcmp(name, "sigma-one")) func = sigma_one;
	if (0 == strcmp(name, "carmichael")) func = carmichael_lambda;
	if (0 == strcmp(name, "mobius")) func = moebius_mu;
	if (0 == strcmp(name, "little-omega")) func = little_omega;
	if (0 == strcmp(name, "big-omega")) func = big_omega;
	if (0 == strcmp(name, "liouville")) func = liouville_lambda;
	if (0 == strcmp(name, "mertens")) func = mertens_m;
	if (0 == strcmp(name, "thue-morse")) func = thue_morse;
	// if (0 == strcmp(name, "")) func =

	if (NULL == func)
	{
		fprintf(stderr, "Oh No, Mr. Bill!  Name %s unknown!\n", name);
		exit(1);
	}

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
	int wrap = 0;  // number of times the phase has wrapped.
	double prev_phase = -1e10;

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
		cpx_exponential_genfunc(val, z, nprec, func);

		// Sum the magnitudes
		cpx_abs(mag, val);
		mpf_add(sum, sum, mag);

		double th = mpf_get_d(theta);
		double mg = mpf_get_d(mag);
		double sm = mpf_get_d(sum);
		sm *= delt;

		// Get the phases, too. atan2 returns a value from -M_PI to +M_PI
		double phase = atan2(cpx_get_im(val), cpx_get_re(val));
		// if (phase < 0.0) phase += 2.0*M_PI;
		if (phase < prev_phase) wrap ++;
		prev_phase = phase;

		printf("%d	%g	%g	%g	%d	%g\n", i, th, mg, sm, wrap, phase);

		mpf_add(theta, theta, delta);
		thet += delt;
		if (2.0*M_PI < thet) break;
	}
}
