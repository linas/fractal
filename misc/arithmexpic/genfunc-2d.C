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

static int divisori(int i) { return divisor(i); }
static double divisor_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, divisori);
	return cabs(g);
}

static double liouv_omega_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, liouville_omega);
	return cabs(g);
}

static double mertens_m_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, mertens_m);
	return cabs(g);
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

// ========================================================
// other stuff.
/* static */ double ploto(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;

      // double t = liouville_lambda (i+1);
      // double t = mangoldt_lambda (i+1);
      // double t = thue_morse (i+1);

	double complex g = ordinary_genfunc(z, totient_phi);
	// double complex g = gpf_normed(z);
	// double complex g = gpf_lambert(z);

	return 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;
}

/* static */ double plot_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z;
	cpx_init(sum);
	cpx_init(z);

	cpx_set_d(z, re_q, im_q);

return 0.0;
// #define PHASE 1
#if PHASE
	// cpx_gpf_ordinary_recip(sum, z, 15);
	cpx_gpf_exponential(sum, z, 20);
	// cpx_gpf_poch_rising(sum, z, 45);
	// cpx_gpf_poch_falling(sum, z, 15);

	double rv = 0.5 + 0.5 * atan2(cpx_get_im(sum), cpx_get_re(sum))/M_PI;
	return rv;
#endif

}

__attribute__((constructor)) void decl_things() {
	DECL_HEIGHT("totient_ord_phase", totient_ord_phase);
	DECL_HEIGHT("totient_exp_phase", totient_exp_phase);
	DECL_HEIGHT("totient_exp_mag", totient_exp_mag);
	DECL_HEIGHT("mobius_exp_mag", mobius_exp_mag);
	DECL_HEIGHT("divisor_exp_mag", divisor_exp_mag);
	DECL_HEIGHT("liouv_omega_exp_mag", liouv_omega_exp_mag);
	DECL_HEIGHT("mertens_m", mertens_m_exp_mag);
	DECL_HEIGHT("thue_morse", thue_morse_exp_mag);
}

// DECL_MAKE_HEIGHT(plot_big);
MAKE_HEIGHT;
