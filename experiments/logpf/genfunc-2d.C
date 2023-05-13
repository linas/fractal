/*
 * Generating functions for miscellaneous logarithmic things
 * 2D phase plot.
 *
 * April 2016; June 2017
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <mp-arith.h>
#include <mp-genfunc.h>
#include <mp-trig.h>

#include <gpf.h>

#include "brat.h"


//  C and C++ is fugnuts insane in complex support.
#define complex _Complex

#define MAX_PREC 1.0e-18
int max_iter = 100000000;

/*
 * Exponential generating function for arithmetic series
 */
double complex exponential_genfunc(long double complex x, unsigned long (*func)(unsigned long))
{
	long double complex sum = 0;
	long double complex xn = x;

	if (cabsl(x) < MAX_PREC) return x;

	for (int n=1; ; n++)
	{
		sum += func(n) * xn;
		xn *= x / ((long double) n+1);
		if (n*cabsl(xn) < MAX_PREC*cabsl(sum)) break;
		if (max_iter < n) break;
	}

	long double scale = expl(-cabsl(x));
	sum *= scale;
	sum *= scale;
// printf("duuude %g sum=%g\n", cabs(x), cabs(sum));

	return sum;
}

/*
 * Exponential generating function for arithmetic series
 */
double complex exp_genfunc_d(long double complex x, double (*func)(long))
{
	long double complex sum = 0;
	long double complex xn = x;
	long double complex maxum = x;
	long double mabs = cabsl(x);

	if (cabsl(x) < MAX_PREC) return x;

	int n;
	for (n=1; ; n++)
	{
		sum += func(n) * xn;
		xn *= x / ((long double) n+1);

		long double sumabs = cabsl(sum);
		if (mabs < sumabs) { mabs = sumabs; maxum = sum; } 
		if (n*cabsl(xn) < MAX_PREC*sumabs) break;
		if (max_iter < n) break;
		// printf("loop at n=%d term=%Lg  sum=%Lg rat=%Lg\n", n, n*cabsl(xn), cabsl(sum), n*cabsl(xn)/cabsl(sum));
	}

// printf("halt at n=%d term=%Lg  sum=%Lg rat=%Lg\n", n, n*cabsl(xn), cabsl(sum), n*cabsl(xn)/cabsl(sum));
// printf(" =======================================\n");

	// long double scale = expl(-cabsl(0.8666666666*x));
	long double scale = expl(-cabsl(x));
	scale *= 1.0e16;
	sum *= scale;
// printf("halt at n=%d sum=%Lg scale=%Lg x=%Lg\n", n, cabsl(sum), scale, cabsl(x));
// printf("duuude %g sum=%g\n", cabs(x), cabs(sum));

	return sum;
}

// ------------------------------------------------------------
//
// Discrete log
// Double precision bombs pretty soon.
//
// double disclog(long n) { return log((double)n+1); }
double disclog(long n) { return log((double)n); }
static double disclog_phase(double re_q, double im_q, int itermax, double param)
{
// max_iter = itermax;
	double complex z = re_q + I * im_q;
	double complex g = exp_genfunc_d(z, disclog);
	double re = creal(g);
	double im = cimag(g);

#if 0
// Draw rulers for the web-page illustration.
double width = 0.03;
double id = im_q/M_PI;
int i = round(id);
if (fabs(id-i) < width) return 1.0;
if (fabs(re_q/M_PI) < width) return 1.0;
#endif
	double rv = 0.5 + 0.5 * atan2(im, re)/M_PI;
	// rv += 0.5;
	// if (1.0 < rv) rv -= 1.0;
	return rv;
}

// Let mpf do the work.
void disc_big(mpf_t ln, long n)
{
	fp_log_ui(ln, n, 165);
}
static double disclog_big_phase(double re_q, double im_q, int itermax, double param)
{
	int nprec = 65;
nprec = itermax;
	mpf_set_default_prec(nprec * 3.322 + 50);

	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	// mpf_t val; mpf_init(val);
	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc_mpf(sum, z, nprec, disc_big);
	// cpx_abs(val, sum);
	// double rv = mpf_get_d(val);

	double re = cpx_get_re(sum);
	double im = cpx_get_im(sum);
	double rv = 0.5 + 0.5 * atan2(im, re)/M_PI;

	cpx_clear(sum);
	cpx_clear(z);
	// mpf_clear(val);
	return rv;
}

static double disclog_mag(double re_q, double im_q, int itermax, double param)
{
	double complex z = re_q + I * im_q;
	double complex g = exp_genfunc_d(z, disclog);
	return cabs(g);
}

// ------------------------------------------------------------

double loggpf(long n) { return log((double)gpf(n)); }
static double loggpf_phase(double re_q, double im_q, int itermax, double param)
{
	double complex z = re_q + I * im_q;
	double complex g = exp_genfunc_d(z, loggpf);
	double re = creal(g);
	double im = cimag(g);
	double rv = 0.5 + 0.5 * atan2(im, re)/M_PI;
	return rv;
}

static double loggpf_mag(double re_q, double im_q, int itermax, double param)
{
	double complex z = re_q + I * im_q;
	double complex g = exp_genfunc_d(z, loggpf);
	return cabs(g);
}

// Let mpf do the work.
void loggpf_big(mpf_t ln, long n)
{
	fp_log_ui(ln, gpf(n), 65);
}
static double loggpf_big_mag(double re_q, double im_q, int itermax, double param)
{
	int nprec = 65;
	mpf_set_default_prec(nprec * 3.322 + 50);

	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);
	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc_mpf(sum, z, nprec, loggpf_big);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

// ------------------------------------------------------------

static double fact_phase(double re_q, double im_q, int itermax, double param)
{
	double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, factor_product);
	double re = creal(g);
	double im = cimag(g);
	double rv = 0.5 + 0.5 * atan2(im, re)/M_PI;
	return rv;
}

static double fact_mag(double re_q, double im_q, int itermax, double param)
{
	double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, factor_product);
	return cabs(g);
}

// ------------------------------------------------------------

static double thue_morse_rev(long n) { return 1.0; }
static double xperiment(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	double complex z = re_q + I * im_q;
	double complex g = exp_genfunc_d(z, thue_morse_rev);
	return cabs(g);
}

// ========================================================


__attribute__((constructor)) void decl_things() {

	DECL_HEIGHT("disclog_phase", disclog_phase);
	DECL_HEIGHT("disclog_big", disclog_big_phase);
	DECL_HEIGHT("disclog_mag", disclog_mag);
	DECL_HEIGHT("loggpf_phase", loggpf_phase);
	DECL_HEIGHT("loggpf_mag", loggpf_mag);
	DECL_HEIGHT("loggpf_big_mag", loggpf_big_mag);

	DECL_HEIGHT("fact_phase", fact_phase);
	DECL_HEIGHT("fact_mag", fact_mag);

	DECL_HEIGHT("xperiment", xperiment);
}

MAKE_HEIGHT;
