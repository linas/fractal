/*
 * Find zeros of function.
 *
 * Perform a survey of a circular domain, serching for zeros
 * in cells of a given angular, radial size.
 *
 * April 2016, October 2016
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <mp-arith.h>
#include <mp-complex.h>
#include <mp-consts.h>
#include <mp-genfunc.h>
#include <mp-misc.h>
#include <mp-trig.h>
#include <mp-zerofind.h>

#include "moebius.h"


// Partition function
void parti_z_mpf(mpf_t res, long n)
{
	mpz_t part; mpz_init(part);
	partition_z(part, n);
	mpf_set_z(res, part);
	mpz_clear(part);
}

void divisor_mpf(mpf_t res, long n)
{
	mpf_set_ui(res, divisor(n));
}

void parti(cpx_t f, cpx_t z, int nprec)
{
	// cpx_exponential_genfunc_mpf(f, z, nprec, parti_z_mpf);
	cpx_exponential_genfunc_mpf(f, z, nprec, divisor_mpf);
}

bool survey_cell(void (*func)(cpx_t f, cpx_t z, int nprec),
                 cpx_t zero,
                 double rguess, double tguess, double cell_size,
                 int ndigits, int nprec)
{
	mp_bitcnt_t bits = ((double) nprec) * 3.3219281 + 50;

	cpx_t a,b,c,d;
	cpx_init2(a, bits);
	cpx_init2(b, bits);
	cpx_init2(c, bits);
	cpx_init2(d, bits);

	cpx_t fa,fb,fc,fd;
	cpx_init2(fa, bits);
	cpx_init2(fb, bits);
	cpx_init2(fc, bits);
	cpx_init2(fd, bits);

	// Set up four corners
	double st = sin(tguess);
	double ct = cos(tguess);
	cpx_set_d(a, rguess*ct, rguess*st);
	double rgd = rguess + cell_size;
	cpx_set_d(b, rgd*ct, rgd*st);
	tguess += cell_size / rguess;
	st = sin(tguess);
	ct = cos(tguess);
	cpx_set_d(c, rgd*ct, rgd*st);
	cpx_set_d(d, rguess*ct, rguess*st);

	// Evaluate at the four corners
	func(fa, a, nprec);
	func(fb, b, nprec);
	func(fc, c, nprec);
	func(fd, d, nprec);

	// Compute contour integral i.e. sum the phases.
	double phase = atan2(cpx_get_im(fa), cpx_get_re(fa));
	if (phase < 0.0) phase += 2.0 * M_PI;
	double sum = phase;
	double hi = phase;
	double lo = phase;

	phase = atan2(cpx_get_im(fb), cpx_get_re(fb));
	if (phase < 0.0) phase += 2.0 * M_PI;
	if (hi < phase) hi = phase;
	if (phase < lo) lo = phase;
	sum += phase;

	phase = atan2(cpx_get_im(fc), cpx_get_re(fc));
	if (phase < 0.0) phase += 2.0 * M_PI;
	if (hi < phase) hi = phase;
	if (phase < lo) lo = phase;
	sum += phase;

	phase = atan2(cpx_get_im(fd), cpx_get_re(fd));
	if (phase < 0.0) phase += 2.0 * M_PI;
	if (hi < phase) hi = phase;
	if (phase < lo) lo = phase;
	sum += phase;

	// We expect the phase to wind around, so that hi and low
	// differ by almost 2pi.
	// We expect the integral to be large, approaching pi,
	// but if we are unlucky, as low as 3pi/8, I guess.
	// if (hi-lo < 4.0 and 0.25*sum < 1.5)
	if (hi-lo < 4.0)
	{
		cpx_clear(a);
		cpx_clear(b);
		cpx_clear(c);
		cpx_clear(d);

		cpx_clear(fa);
		cpx_clear(fb);
		cpx_clear(fc);
		cpx_clear(fd);

		return false;
	}

	// Now, a will be the center of the rectangle.
	cpx_add(a, a, b);
	cpx_add(a, a, c);
	cpx_add(a, a, d);
	cpx_div_ui(a, a, 4);

	cpx_t e1, e2;
	cpx_init2(e1, bits);
	cpx_init2(e2, bits);
	cpx_sub(e1, b, a);
	cpx_sub(e2, d, a);

	int rc = cpx_find_zero(zero, func, a, e1, e2, ndigits, nprec);

if (rc) {printf("duuude found noothing\n");}

	// if rc is not zero, then nothing was found
	if (rc)
	{
		cpx_clear(a);
		cpx_clear(b);
		cpx_clear(c);
		cpx_clear(d);

		cpx_clear(fa);
		cpx_clear(fb);
		cpx_clear(fc);
		cpx_clear(fd);

		cpx_clear(e1);
		cpx_clear(e2);
		return false;
	}

	printf("---------------\n");
	cpx_prt("zero = ", zero); printf("\n");

	mpf_t r, t, pi;
	mpf_init2(r, bits);
	mpf_init2(t, bits);
	mpf_init2(pi, bits);

	cpx_abs(r, zero);
	fp_prt("r = ", r); printf("\n");

	fp_arctan2(t, zero->im, zero->re, nprec);
	fp_pi(pi, nprec);
	mpf_div(t, t, pi);

	fp_prt("theta/pi = ", t); printf("\n");

#define CHECK_RESULT 1
#ifdef CHECK_RESULT
	cpx_t check;
	cpx_init(check);
	func(check, zero, nprec);
	double eps_r = cpx_get_re(check);
	double eps_i = cpx_get_im(check);
	double eps = sqrt(eps_r*eps_r + eps_i*eps_i);
	printf("fun at zero = %g\n", eps);
	cpx_clear(check);
#endif

	mpf_clear(r);
	mpf_clear(t);
	mpf_clear(pi);

	printf("\n");
	fflush(stdout);

	cpx_clear(a);
	cpx_clear(b);
	cpx_clear(c);
	cpx_clear(d);

	cpx_clear(fa);
	cpx_clear(fb);
	cpx_clear(fc);
	cpx_clear(fd);

	cpx_clear(e1);
	cpx_clear(e2);
	return true;
}

void survey(void (*func)(cpx_t f, cpx_t z, int nprec),
            double rmax, double cell_size,
            int ndigits, int nprec)
{
	mp_bitcnt_t bits = ((double) nprec) * 3.3219281 + 50;
	cpx_t zero;
	cpx_init2(zero, bits);
	for (double r=1.0; r<rmax; r += cell_size)
	{
		double step = cell_size / r;
		for (double t=0.0; t < M_PI; t += step)
		{
			survey_cell(func, zero, r, t, cell_size, ndigits, nprec);
		}
	}
	cpx_clear(zero);
}

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s <rmax> <cell-size>\n", argv[0]);
		exit(1);
	}

	int nprec = 80;
	mp_bitcnt_t bits = ((double) nprec) * 3.3219281 + 50;
	mpf_set_default_prec(bits);

	double rmax = atof(argv[1]);
	double cell_size = atof(argv[2]);
	survey(parti, rmax, cell_size, 40, nprec);

	printf("Done!\n");
}
