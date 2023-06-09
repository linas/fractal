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
		xn *= x / ((long double) n+1);
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

/* static */ double ploto(double re_q, double im_q, int itermax, double param)
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

#define HYPERBOLIC
#ifdef HYPERBOLIC

// #define Q_DISK
#ifdef Q_DISK
	// This code converts from the q-disk (where q == nome) to upper half-plane.
	// That is, the input coords are the nome, the output are half-plane coords.
	double mod = sqrt(re_q * re_q + im_q * im_q);
	if (1.0 <= mod) return 1.0;
	double hpy = -log(mod) / M_PI;
	double hpx = atan2(im_q, re_q) / M_PI;
#endif // Q_DISK

#define POINCARE_DISK
#ifdef POINCARE_DISK
	// This code converts from poincare disk to upper half-plane
	// Let pdx, pdy be coordinates on the Poincare disk
	// Rotate by i
	double pdx = -im_q;
	double pdy = re_q;
	double rsq = pdx*pdx+pdy*pdy;
	if (1.0 <= rsq) return 1.0;

	// Let hpx, hpy be coords of point in the Poincare upper half-plane ...
	double hpd = pdx*pdx + (1.0-pdy)*(1.0-pdy);
	double hpx = 2.0 * pdx / hpd;
	double hpy = (1.0 - pdx*pdx - pdy*pdy) / hpd;
#endif // POINCARE_DISK

// #define POLAR_PARAMETRIC
#ifdef POLAR_PARAMETRIC
	// As hpx runs from -inf to +inf so theta runs from 0 to 2pi
	double theta = atan2 (hpx, 1.0) + M_PI;

	// The parametric curve is const = sqrt(r) sin(0.5 theta)
	// For const == 1, this encloses the primary zero-free-lane.
	// Ergo, we are interested in the map where const == 1/hpy
	double par = sin(0.5 * theta);
	par = par*par;
	double are = par / (hpy*hpy);
#endif // POLAR_PARAMETRIC

#define HOKEY
#ifdef HOKEY
	// Perform a hokey mapping of upper-half-plane to entire complex plane
	// Compared to the others above, this one looks the most like the
	// original deal (when used with Poincare disk map.)
	double are = 1.0 / (hpy*hpy);
	// are *= are*are*are;

	if (hpy<=0.0) hpy = 1.0e-16;
	are = exp(4.0/hpy) - 1.0;

	double theta = M_PI * tanh(hpx);
#endif

	// Avoid long-running calculations.
	if (itermax < are) return 0.5;

	// Finally, cartesian coordinates on the complex z-plane
	double czx = are * cos(theta);
	double czy = are * sin(theta);

	re_q = czx;
	im_q = czy;

#endif // HYPERBOLIC

// #define HYPERBOLOID
#ifdef HYPERBOLOID
	// Let pdx, pdy be coordinates on the Poincare disk
	double pdx = re_q;
	double pdy = im_q;
	double rsq = pdx*pdx+pdy*pdy;
	if (1.0 <= rsq) return 1.0;

	// Let lbx, bly be coordinates on the hyperboloid
	// double tee = (1.0 + rsq) / (1.0 - rsq);
	// double scale = 50.0;
	double scale = 1.0;
	double blx = scale * 2.0 * pdx / (1.0 - rsq);
	double bly = scale * 2.0 * pdy / (1.0 - rsq);

	re_q = blx;
	im_q = bly;

	double br = sqrt(blx*blx+bly*bly);
	if (0.0 < br)
	{
		double bsin = bly / br;
		double bcos = blx / br;
		re_q = (exp(4.0*br) - 1.0)*bcos;
		im_q = (exp(4.0*br) - 1.0)*bsin;
	}
	else
	{
		re_q = im_q = 0.0;
	}
// printf("duuude blx = %g %g re= %g %g\n", blx, bly, re_q, im_q);
	double are = sqrt(re_q*re_q + im_q*im_q);
	if (itermax < are) return 0.5;
#endif

	// ----------------------

	cpx_t sum, z;
	cpx_init(sum);
	cpx_init(z);
	cpx_set_d(z, re_q, im_q);

// #define ODF_PHASE 1
#if ODF_PHASE
	// cpx_gpf_ordinary_shift(sum, z, itermax, 12);
	cpx_gpf_ordinary_derivative(sum, z, itermax, 12);

	double rv = 0.5 + 0.5 * atan2(cpx_get_im(sum), cpx_get_re(sum))/M_PI;
	return rv;
#endif

// #define PHASE 1
#if PHASE
	// cpx_gpf_ordinary_recip(sum, z, 15);
	cpx_gpf_exponential(sum, z, 20);
	// cpx_gpf_poch_rising(sum, z, 45);
	// cpx_gpf_poch_falling(sum, z, 15);

	double rv = 0.5 + 0.5 * atan2(cpx_get_im(sum), cpx_get_re(sum))/M_PI;
	return rv;
#endif

// #define EXPO 1
#if EXPO
	cpx_gpf_exponential(sum, z, 20);
	// cpx_gpf_sine(sum, z, 20);
	// cpx_gpf_exponential_shift(sum, z, itermax, 25);
	// cpx_gpf_exponential_newton(sum, z, itermax, 25);

	// extract
	mpf_t val;
	mpf_init(val);
	cpx_abs(val, sum);

	double rv = mpf_get_d(val);
// rv = cpx_get_re(sum);

	// Divide by z for plotting.
	double r = sqrt(re_q*re_q + im_q*im_q);
	double lr = log(r);
	rv *= lr / r;
	// double regul = exp(-1.0/(r*r));
	// rv *= (1.0-regul) + regul * lr*lr / (r*r);

	return rv;
#endif

// #define DIFEXPO 1
#if DIFEXPO
	// This is mostly junk, it takes the difference between the
	// expo and its square, and its .. pretty but its junk.
	complex double zee = re_q + I*im_q;
	complex double zp = zee*zee;
	double pre = creal(zp);
	double pim = cimag(zp);
	double are = sqrt(pre*pre + pim*pim);
	if (itermax < are) return 0.5;

	cpx_set_d(z, re_q, im_q);
	cpx_gpf_exponential(sum, z, 40);

	cpx_t sumzsq;
	cpx_init(sumzsq);
	cpx_set_d(z, pre, pim);
	cpx_gpf_exponential(sumzsq, z, 40);

	#ifdef SQUARE
		cpx_t sumsq;
		cpx_init(sumsq);
		cpx_mul(sumsq, sum, sum);
		cpx_sub(sum, sumzsq, sumsq);
	#endif
	cpx_div_ui(sumzsq, sumzsq, 2);
	cpx_sub(sum, sum, sumzsq);

	// extract
	mpf_t val;
	mpf_init(val);
	cpx_abs(val, sum);

	double rv = mpf_get_d(val);

	// Divide by z for plotting.
	double r = sqrt(re_q*re_q + im_q*im_q);
// r = r*r;
	double lr = log(r);
	rv *= lr / r;

	return rv;
#endif
// #define RANDY 1
#if RANDY
	cpx_random_exponential_shift(sum, z, itermax, 25);

	// extract
	mpf_t val;
	mpf_init(val);
	cpx_abs(val, sum);

	double rv = mpf_get_d(val);

	// Divide by z for plotting.
	double r = sqrt(re_q*re_q + im_q*im_q);
	double lr = log(r);
	rv *= (lr*lr) / r;

	return rv;
#endif

#define RECIP 1
#ifdef RECIP

	// #define UN_CIRCLE 1
	#ifdef UN_CIRCLE
		// printf("duuude in= %f %f \n", re_q, im_q);
		double theta = M_PI * im_q;

		#ifdef INSCRIBE
			double x = re_q;
			double y = sin(0.5*theta);
			// if (0.9*param < x*x*y and x*x*y < 1.1*param) return 0.0;
			x = pow(x, 1.5);
			if (0.9*param < x*y and x*y < 1.1*param) return 0.0;
		#endif


		#if LINEAR_CENTER_LINE
			double rr = itermax;
			rr = exp(rr * M_LN2);  // pow (2, itermax * re_q)
			rr += param * re_q; // left-right offsets.
		#endif
		double rr = itermax + param * re_q;
		rr = exp(rr * M_LN2);  // pow (2, param * re_q)

		im_q = rr*sin (theta);
		re_q = rr*cos (theta);
		cpx_set_d(z, re_q, im_q);

		// printf("duuude              out= %f %f \n", re_q, im_q);
	#endif // UN_CIRCLE

	cpx_gpf_exponential_recip(sum, z, 25);
	// extract
	mpf_t val;
	mpf_init(val);
	cpx_abs(val, sum);

	double rv = mpf_get_d(val);

	// Divide by z for plotting.
	double r = sqrt(re_q*re_q + im_q*im_q);
	double lr = log(r+1);
	double llr = log(lr+1);
	// rv *= lr * lr * lr;
	// rv /= (lr * lr * lr);  // used for the off-by-one pictures
	// rv *= llr * llr * llr;
	// rv *= lr * lr * lr *lr;
	// rv *= lr * lr * lr *lr *lr;  // used for the main pictures
	rv *= llr * llr * llr *llr *lr;

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

/* static */ double plot_diri(double re_q, double im_q, int itermax, double param)
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
