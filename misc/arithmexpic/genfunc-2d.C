/*
 * Generating functions for miscellaneous arithmetic series
 * 2D phase plot.
 *
 * April 2016
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include <mp-zeta.h>

#include <gpf.h>
#include <modular.h>
#include <moebius.h>
#include <totient.h>

#include "brat.h"
#include "genfunc.h"



//  C and C++ is fugnuts insane in complex support.
#define complex _Complex

#define MAX_PREC 1.0e-18
int max_iter = 100000000;

/*
 * Ordinary generating function for arithmetic series
 */
double complex ordinary_genfunc(double complex x, int (*func)(int))
{
	double complex sum = 0;
	double complex xn = x;

	if (cabs(x) < 1.0e-16) return x;
	if (0.9999999 <= cabs(x)) return 0.0;

	for (int n=1; ; n++)
	{
		sum += func(n) * xn;
		xn *= x;
		if (n*cabs(xn) < MAX_PREC*cabs(sum)) break;
		if (max_iter < n) break;
	}

	return sum;
}

/*
 * Exponential generating function for arithmetic series
 */
double complex exponential_genfunc(long double complex x, int (*func)(int))
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
// printf("duuude %g sum=%g\n", cabs(x), cabs(sum));

	return sum;
}

/*
 * Exponential generating function for arithmetic series
 */
double complex exp_genfunc_d(long double complex x, double (*func)(int))
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
// printf("duuude %g sum=%g\n", cabs(x), cabs(sum));

	return sum;
}

double complex lambert_genfunc(double complex x, int (*func)(int))
{
	double complex sum = 0;
	double complex xn = x;

	if (cabs(x) < 1.0e-16) return x;
	if (0.9999999 <= cabs(x)) return 0.0;

	for (int n=1; ; n++)
	{
		sum += func(n) * xn / (1.0 - xn);
		xn *= x;
		if (n*cabs(xn) < MAX_PREC*cabs(sum)) break;
		if (max_iter < n) break;
	}

	return sum;
}

// ======================================================================
// Totient stuff...

static double totient_ord_phase(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = ordinary_genfunc(z, totient_phi);
	return 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;
}

static double totient_exp_phase(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, totient_phi);
	return 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;
}

static double totient_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, totient_phi);
	return cabs(g);
}

static double mobius_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, moebius_mu);
	return cabs(g);
}

static double mobius_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);

	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc(sum, z, 25, moebius_mu);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	return rv;
}

static int divisori(int i) { return divisor(i); }
static double divisor_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, divisori);
	return cabs(g);
}

static double divisor_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);

	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc(sum, z, 25, divisori);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	return rv;
}

static double liouv_omega_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, liouville_omega);
	return cabs(g);
}

static double liouv_lambda(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, liouville_lambda);
	return cabs(g);
}

static double mertens_m_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, mertens_m);
	return cabs(g);
}

double mango(int n) { return mangoldt_lambda_cached(n); }
static double mangoldt_lambda_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = exp_genfunc_d(z, mango);
	return cabs(g);
}

static double exp_mangoldt_lambda_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, exp_mangoldt_lambda);

	double rv = cabs(g);
	double r = cabs(z);
	double lr = log(r+1.0);
	rv /= lr*lr;

	return rv;
}

static double exp_mango_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);

	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc(sum, z, 25, exp_mangoldt_lambda);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	double r = sqrt(re_q*re_q + im_q*im_q);
	double lr = log(r+1.0);
	rv /= lr*lr;
	return rv;
}

static int thue_morse(int n)
{
   if (0 == n) return 0;
   if (1 == n) return 1;
   if (0 == n%2) return thue_morse (n/2);
   return (1-thue_morse ((n-1)/2));
}

static double thue_morse_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, thue_morse);
	return cabs(g);
}

static double thue_morse_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);

	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc(sum, z, 25, thue_morse);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	return rv;
}

// static double thue_morse_recip(int n) { return 1.0 / (1.0 + thue_morse(n)); }
static double thue_morse_rev(int n) { return 1.0 - thue_morse(n); }
static double xperiment(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = exp_genfunc_d(z, thue_morse_rev);
	return cabs(g);
}

// ========================================================
// other stuff.
/* static */ double ploto(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;

	double complex g = ordinary_genfunc(z, totient_phi);
	// double complex g = gpf_normed(z);
	// double complex g = gpf_lambert(z);

	return 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;
}


__attribute__((constructor)) void decl_things() {
	DECL_HEIGHT("totient_ord_phase", totient_ord_phase);
	DECL_HEIGHT("totient_exp_phase", totient_exp_phase);
	DECL_HEIGHT("totient_exp_mag", totient_exp_mag);
	DECL_HEIGHT("mobius_exp_mag", mobius_exp_mag);
	DECL_HEIGHT("mobius_big", mobius_big);
	DECL_HEIGHT("divisor_exp_mag", divisor_exp_mag);
	DECL_HEIGHT("divisor_big", divisor_big);
	DECL_HEIGHT("liouv_omega_exp_mag", liouv_omega_exp_mag);
	DECL_HEIGHT("liouv_lambda", liouv_lambda);
	DECL_HEIGHT("mertens_m", mertens_m_exp_mag);
	DECL_HEIGHT("mangoldt_lambda", mangoldt_lambda_exp_mag);
	DECL_HEIGHT("exp_mangoldt_lambda", exp_mangoldt_lambda_exp_mag);
	DECL_HEIGHT("exp_mango_big", exp_mango_big);
	DECL_HEIGHT("thue_morse", thue_morse_exp_mag);
	DECL_HEIGHT("thue_morse_big", thue_morse_big);
	DECL_HEIGHT("xperiment", xperiment);
}

// DECL_MAKE_HEIGHT(plot_big);
MAKE_HEIGHT;
