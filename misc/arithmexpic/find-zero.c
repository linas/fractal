/*
 * Find zeros of function
 *
 * April 2016, October 2016
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <mp-arith.h>
#include <mp-zerofind.h>
#include "genfunc.h"

// Partition function
void parti_z_mpf(mpf_t res, long n)
{
	mpz_t part; mpz_init(part);
	partition_z(part, n);
	mpf_set_z(res, part);
	mpz_clear(part);
}

void parti(cpx_t f, cpx_t z, int nprec)
{
	cpx_exponential_genfunc_mpf(f, z, nprec, parti_z_mpf);
}

bool survey_cell(void (*func)(cpx_t f, cpx_t z, int nprec),
                 double rguess, double tguess, double cell_size, int nprec)
{
	mp_bitcnt_t bits = 10 + ((double) prec) * 3.3219281;

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
	tguess *= 2.0 * M_PI;
	double st = sin(tguess);
	double ct = cos(tguess);
	cpx_set_d(a, rguess*ct, rguess*st);
	double rgd = rguess + cell_size;
	cpx_set_d(b, rgt*ct, rgt*st);
	tguess += 2.0 * M_PI * cell_size / rguess
	st = sin(tguess);
	ct = cos(tguess);
	cpx_set_d(c, rgt*ct, rgt*st);
	cpx_set_d(d, rguess*ct, rguess*st);

	// Evaluate at the four corners
	func(fa, a, nprec);
	func(fb, b, nprec);
	func(fc, c, nprec);
	func(fd, d, nprec);

	// Compute contour integral i.e. sum the phases.
	double phase = atan2(cimag(fa), creal(fa));
	double sum = phase;
	double hi = phase;
	double lo = phase;
	phase = atan2(cimag(fb), creal(fb));
	sum += phase;
	if (hi < phase) hi = phase;
	if (phase < lo) lo = phase;

	phase = atan2(cimag(fc), creal(fc));
	sum += phase;
	if (hi < phase) hi = phase;
	if (phase < lo) lo = phase;

	phase = atan2(cimag(fd), creal(fd));
	sum += phase;
	if (hi < phase) hi = phase;
	if (phase < lo) lo = phase;

	// We expect the phase to wind around, so that hi and low
	// differ by almost 2pi.
	if (hi-lo < 3.0) return false;

	// We expect the integral to be large, approaching 2pi.
	if (phase < 3.0) return false;

	cpx_t e1, e2, zero;
	cpx_init2(e1, bits);
	cpx_init2(e2, bits);
	cpx_init2(zero, bits);
	cpx_minus(e1, b, a);
	cpx_minus(e2, c, a);

	int rc = cpx_find_zero(zero, func, a, e1, e2, 20, nprec);

	// if rc is not zero, then nothing was found
	if (rc) return;

	double re = cpx_get_re(zero);
	double im = cpx_get_im(zero);

	double r = sqrt(re*re + im*im);
	double t = atan2(im, re) / M_PI;

	printf("z = %-16.14g * exp(i pi %-16.14g)", r, t);
	printf(" = %-16.14g + I %-16.14g\n", re, im);
	fflush(stdout);

#define CHECK_RESULT 1
#ifdef CHECK_RESULT
	cpx_t check;
	cpx_init(check);
	func(check, zero, 20);
	double eps_r = cpx_get_re(check);
	double eps_i = cpx_get_im(check);
	double eps = sqrt(eps_r*eps_r + eps_i*eps_i);
	printf("fun at zero = %g\n", eps);
#endif
}

void survey(void (*func)(cpx_t f, cpx_t z, int nprec),
            double rmax, double cell_size, int nprec)
{
	for (double r=cell_size; r<rmax; r += cell_size)
	{
		for (double t=0.0; t < 0.5; t += cell_size/r)
		{
			survey_cell(func, r, t, cell_size, nprec);
			restart = false;
		}
	}
}

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s <rmax> <cell-size>\n", argv[0]);
		exit(1);
	}

	double rmax = atof(argv[1]);
	double cell_size = atof(argv[2]);
	survey(parti, rmax, cell_size, 60);
}
