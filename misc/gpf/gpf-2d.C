/*
 * Generating functions for greatest prime factors.
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
#include "brat.h"
#include "gpf-gen-bignum.h"

//  C and C++ is fugnuts insane in complex support.
#define complex _Complex

#define MAX_PREC 1.0e-18
int max_iter = 100000000;

/*
 * Ordinary generating function for the greatest common factor.
 */
double complex gpf_ordinary(double complex x)
{
	double complex sum = 0;
	double complex xn = x;

	if (cabs(x) < 1.0e-16) return x;
	if (0.9999999 <= cabs(x)) return 0.0;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn;
		xn *= x;
		if (n*cabs(xn) < MAX_PREC*cabs(sum)) break;
		if (max_iter < n) break;
	}

	return sum;
}

double complex gpf_normed(double complex x)
{
	double complex sum = 0;
	double complex xn = x;

	if (cabs(x) < 1.0e-16) return x;
	if (0.9999999 <= cabs(x)) return 0.0;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn / ((double) n);
		xn *= x;
		if (cabs(xn) < MAX_PREC*cabs(sum)) break;
		if (max_iter < n) break;
	}

	return sum;
}

/*
 * Exponential generating function for the greatest common factor.
 */
double complex gpf_exponential(long double complex x)
{
	long double complex sum = 0;
	long double complex xn = x;

	if (cabsl(x) < MAX_PREC) return x;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn;
		xn *= x / ((long double) n);
		if (n*cabsl(xn) < MAX_PREC*cabsl(sum)) break;
		if (max_iter < n) break;
	}

	long double scale = expl(-cabsl(x));
	sum *= scale;
// printf("duuude %g sum=%g\n", cabs(x), cabs(sum));

	return sum;
}

double complex gpf_lambert(double complex x)
{
	double complex sum = 0;
	double complex xn = x;

	if (cabs(x) < 1.0e-16) return x;
	if (0.9999999 <= cabs(x)) return 0.0;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * xn / (1.0 - xn);
		xn *= x;
		if (n*cabs(xn) < MAX_PREC*cabs(sum)) break;
		if (max_iter < n) break;
	}

	return sum;
}

static double ploto(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
   double complex z = re_q + I * im_q;

	// double complex g = gpf_ordinary(z);
	double complex g = gpf_exponential(z);
	// g *= cexp(-z);
	// g *= exp(-cabs(z)) / cabs(z);
	// double complex g = gpf_normed(z);
	// double complex g = gpf_lambert(z);

	// double rv = cabs(g);
	// double r = sqrt(re_q*re_q + im_q*im_q);
	// rv /= r;
	// rv /= r*r/log(r);
	// return rv;

	// return cabs(g);
	// return creal(g);

	return 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;
}

static double plot_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z;
	cpx_init(sum);
	cpx_init(z);

	cpx_set_d(z, re_q, im_q);

// #define PHASE 1
#if PHASE
	// cpx_gpf_ordinary_recip(sum, z, 15);
	// cpx_gpf_exponential(sum, z, 20);
	cpx_gpf_poch_rising(sum, z, 45);
	// cpx_gpf_poch_falling(sum, z, 15);

	double rv = 0.5 + 0.5 * atan2(cpx_get_im(sum), cpx_get_re(sum))/M_PI;
	return rv;
#endif

#define EXPO 1
#if EXPO
	// cpx_gpf_exponential(sum, z, 20);
	cpx_gpf_exponential_d(sum, z, itermax, 25);

	// extract
	mpf_t val;
	mpf_init(val);
	cpx_abs(val, sum);

	double rv = mpf_get_d(val);
// rv = cpx_get_re(sum);

	// Divide by z for plotting.
	double r = sqrt(re_q*re_q + im_q*im_q);
	double lr = log(r);
	// rv /= r*r / (lr*lr);
	rv *= (1.0 - exp(-1.0/(r*r))) * lr*lr / (r*r);

	return rv;
#endif

// #define RECIP 1
#ifdef RECIP
	cpx_gpf_exponential_recip(sum, z, 15);
	// extract
	mpf_t val;
	mpf_init(val);
	cpx_abs(val, sum);

	double rv = mpf_get_d(val);

	// Divide by z for plotting.
	double r = sqrt(re_q*re_q + im_q*im_q);
	double lr = log(r);
	rv /= (lr * lr * lr);

	return rv;
#endif

// #define ESS 1
#ifdef ESS
	cpx_t s;
	cpx_init(s);

	cpx_set_d(s, param, 0.0);

	cpx_gpf_exponential_s(sum, z, s, 15);
	// extract
	mpf_t val;
	mpf_init(val);
	cpx_abs(val, sum);

	double rv = mpf_get_d(val);

rv = cpx_get_re(sum);
	// Divide by z for plotting.
	double r = sqrt(re_q*re_q + im_q*im_q);
	double lr = log(r);
	// rv /= (lr * lr);

	// Standard S=1 normalization.
	rv /= r*r / (lr*lr);

	return rv;
#endif

// #define POCH 1
#ifdef POCH
	cpx_t fal;
	cpx_init(fal);
	cpx_gpf_poch_rising(sum, z, 45);
	// cpx_gpf_poch_falling(sum, z, 35);
	// cpx_gpf_poch_falling(fal, z, 25);
	// cpx_sub(sum, sum, fal);

	// extract
	mpf_t val;
	mpf_init(val);
	cpx_abs(val, sum);

	double rv = mpf_get_d(val);
	double r = sqrt(re_q*re_q + im_q*im_q);
	// double lr = log(r);
// double gre = cpx_get_re(sum);
// double gim = cpx_get_im(sum);
// printf("duude re=%g im=%g gre=%g gim=%g\n", re_q, im_q, gre, gim);

	// rv *= exp(-2.0*sqrt(r));
	rv /= r;
	// rv *= 5.0;
#if 0
double lv = log(rv);
if (lv < 0.0) lv = 0.0;
printf("duude rv=%g scale=%g\n", rv, lv/lr);
#endif
// printf("duude r=%g rv=%g \n", r, rv);

	return rv;
#endif
}

static double plot_diri(double re_q, double im_q, int itermax, double param)
{
static int cnt=0;
int id = ++cnt;
	// discard outside of the unit circle.
	if (1.0 <= re_q*re_q + im_q*im_q) return 0.0;

	// Map the inside of a unit circle to the right-hand complex
	// half-plane.

	// if q in circle, and z in upper half plane, then
	// q = (z-1)/(z+1)  or z = (1+q) / (1-q)
	double complex q = re_q + I * im_q;
	double complex z = (1.0 + q) / (1.0 - q);

	// Next, we want to rotate by 90 and offset.
	// Thus, z becomes s so that |s| > offset
	z *= -I;
	z += 4.0;

	// Finally, avoid travelling too far up the imaginary axis,
	// as this hinders convergence.
	if (4.0 < fabs(cimag(z))) return 0.0;

	if (3.6 >= creal(z)) return 0.0;

	cpx_t sum, ess;
	cpx_init(sum);
	cpx_init(ess);

time_t start = time(NULL);
printf("Start pix=%d start work on %g %g\n", id, creal(z), cimag(z));
	cpx_set_d(ess, creal(z), cimag(z));

	cpx_gpf_dirichlet(sum, ess, 15);
	// cpx_borwein_zeta(sum, ess, 15);
time_t stop = time(NULL);

	double rv = 0.5 + 0.5 * atan2(cpx_get_im(sum), cpx_get_re(sum))/M_PI;
printf("Done pix=%d done took %lu on %g %g val=%g\n", id, stop-start, creal(z), cimag(z), rv);
	return rv;
}

// DECL_MAKE_HEIGHT(ploto);
DECL_MAKE_HEIGHT(plot_big);
// DECL_MAKE_HEIGHT(plot_diri);
