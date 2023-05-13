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

#include <mp-arith.h>
#include <mp-genfunc.h>
#include <mp-trig.h>

#include <gpf.h>
#include <isqrt.h>
#include <moebius.h>
#include <necklace.h>
#include <totient.h>

#include "brat.h"



//  C and C++ is fugnuts insane in complex support.
#define complex _Complex

#define MAX_PREC 1.0e-18
int max_iter = 100000000;

/*
 * Ordinary generating function for arithmetic series
 */
double complex ordinary_genfunc(double complex x, long (*func)(long))
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
double complex exponential_genfunc(long double complex x, long (*func)(long))
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
double complex exp_genfunc_d(long double complex x, double (*func)(long))
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

double complex lambert_genfunc(double complex x, long (*func)(long))
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
	long double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, totient_phi);
	return 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;
}

static double totient_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	long double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, totient_phi);
	return cabs(g);
}

static double totient_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);

	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc(sum, z, 25, totient_phi);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	double r = sqrt(re_q*re_q + im_q*im_q);
	double lr = log(r+1.0);
	rv /= lr;

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

static double carmichael(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	long double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, carmichael_lambda);
	return cabs(g);
}

static double carmichael_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);

	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc(sum, z, 25, carmichael_lambda);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	double r = sqrt(re_q*re_q + im_q*im_q);
	double lr = log(r+1.0);
	rv /= lr *lr;

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

static double mobius_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
// #define UNCIRCLE
#ifdef UNCIRCLE
	double tmp = re_q;
	re_q = 1.0 - im_q;
	im_q = tmp;

	// max_iter = itermax;
	max_iter = 100000;
	double theta = M_PI * im_q;
	double rr = itermax + param * re_q;
	rr = exp(rr * M_LN2);  // pow (2, itermax + param * re_q)
// printf("duuude re=%g im=%g r = %g\n", re_q, im_q, rr);
	im_q = rr*sin (theta);
	re_q = rr*cos (theta);
#endif

	long double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, moebius_mu);
	return 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;
	// return cabs(g);
}

static double mobius_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);

	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc(sum, z, 25, moebius_mu);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

static double divisor_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	long double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, divisor);
	return cabs(g);
}

static double divisor_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);

	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc(sum, z, 25, divisor);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

/* Unwind the circle into a rectangle */
static double divisor_uncircle(double re_q, double im_q, int itermax, double param)
{
	double tmp = re_q;
	re_q = 1.0 - im_q;
	im_q = tmp;

	// max_iter = itermax;
	max_iter = 100000;
	double theta = M_PI * im_q;
	double rr = itermax + param * re_q;
	rr = exp(rr * M_LN2);  // pow (2, itermax + param * re_q)
// printf("duuude re=%g im=%g r = %g\n", re_q, im_q, rr);
	im_q = rr*sin (theta);
	re_q = rr*cos (theta);

	long double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, divisor);

#ifdef RULER
double dr = 0.1*(log(rr) - 2);
int ir = rr;
ir = ir - 10*(ir%10);
if (rr-dr < ir and ir < rr+dr) return 0.0;
#endif

	return 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;
	// return cabs(g);
// #define STYLIZE
#ifdef STYLIZE
	// phi the phase runs from zero to 1.
	double phi = 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;

	double yellow = 0.75;
	double ywidth = 0.03;
	double gwidth = 0.03;
	double side = 0.0;
	if (1.0-ywidth < phi) side = yellow;
	if (0.5-gwidth < phi and phi < 0.5+gwidth) side=yellow;
	if (phi< ywidth) side = yellow;

	return side;
#endif

}

/* The below attempts to remap a pie slice of the circle into
 * the whole circle, as suggested by the dyadic self-similarity
 * However, it is a total failure as such. The mapping is not
 * holomorphic, its not even close. It sucks.
 */
static double divisor_twist(double re_q, double im_q, int itermax, double param)
{
	int prec = 35;
	mp_bitcnt_t bits = ((double) prec) * 3.322 + 50;

	cpx_t sum, z; cpx_init2(sum, bits); cpx_init2(z, bits);
	mpf_t val; mpf_init2(val, bits);

	cpx_set_d(z, re_q, im_q);

	cpx_exponential_twist(sum, z, prec, divisor);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

static long sigma1(long i) { return sigma(i,1); }
static double sigma_one(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	long double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, sigma1);

	double rv = cabs(g);
	double r = sqrt(re_q*re_q + im_q*im_q);
	double lr = log(r+1.0);
	rv /= lr;
	return rv;
}

static long sigma2(long i) { return sigma(i,2); }
static double sigma_two(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	long double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, sigma2);
	double rv = cabs(g);
	double r = sqrt(re_q*re_q + im_q*im_q);
	double lr = log(r+1.0);
	rv /= lr;
	rv /= lr;
	return rv;
}

static double little_omega_big(double re_q, double im_q, int itermax, double param)
{
	int nprec = 85;
	mpf_set_default_prec(nprec * 3.322 + 50);

	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);
	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc(sum, z, 25, little_omega);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

static long unitary_divisor(long n) { return 1<<little_omega(n); }
static double unitary_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);
	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc(sum, z, 25, unitary_divisor);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

static double big_omega_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	long double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, big_omega);
	return cabs(g);
}

static double big_omega_big(double re_q, double im_q, int itermax, double param)
{
	int nprec = 45;
	mpf_set_default_prec(nprec * 3.322 + 50);

	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);
	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc(sum, z, 25, big_omega);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

static double liouv_lambda(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	long double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, liouville_lambda);
	return cabs(g);
}

static double mertens_m_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	long double complex z = re_q + I * im_q;
	double complex g = exponential_genfunc(z, mertens_m);
	return cabs(g);
}

// ----------------------------------------------------------------
// von Mangoldt fuction
double mango(long n) { return mangoldt_lambda_cached(n); }
static double mangoldt_lambda_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	double complex z = re_q + I * im_q;
	double complex g = exp_genfunc_d(z, mango);
	return cabs(g);
}

// Well we need the logarithm; let mpf do the work.
void mango_big(mpf_t ln, long n)
{
	fp_log_ui(ln, exp_mangoldt_lambda(n), 25);
}
static double mangoldt_lambda_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);
	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc_mpf(sum, z, 25, mango_big);
	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

static double exp_mangoldt_lambda_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	long double complex z = re_q + I * im_q;
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

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

static double thue_morse_exp_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	long double complex z = re_q + I * im_q;
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

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

static double isqrt_big(double re_q, double im_q, int itermax, double param)
{
	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);

	cpx_set_d(z, re_q, im_q);

	cpx_exponential_genfunc(sum, z, 25, isqrt);

#if 0
	mpf_t gabs; mpf_init(gabs);
	cpx_abs(gabs, z);
	mpf_sqrt(gabs, gabs);
	cpx_div_mpf(sum, sum, gabs);
#endif

	cpx_abs(val, sum);
	double rv = mpf_get_d(val);

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

long foop(long n)
{
	static long maxn=0, maxv;
	long rv = partition(n);
if (rv == 0) return 1;
	if (maxn < n) {
		maxn = n;
		maxv = rv;
		int b = 0;
		long ss = rv;
		while (0 < ss) { b++; ss>>=1; }
		printf("partition n=%ld v=%ld bits=%d\n", n, maxv, b);
	}
	return rv;
}

void parti_mpf(mpf_t res, long n)
{
	unsigned __int128 rv = partitionll(n);

	static bool overflow = false;
	if (0 == rv and not overflow)
	{
		overflow = true;
		printf("Partition Overflow!!!!!!\n");
	}
	if (rv == 0) { mpf_set_ui(res, 1); return; }

	unsigned __int128 mask = 1UL<<30;
	mask <<= 34;
	mask -= 1;
	unsigned long lo = rv & mask;
	unsigned long hi = rv >> 64;
	mpf_set_ui(res, hi);
	mpf_mul_ui(res, res, 1UL<<30);
	mpf_mul_ui(res, res, 1UL<<34);
	mpf_add_ui(res, res, lo);
}

void parti_z_mpf(mpf_t res, long n)
{
	mpz_t part; mpz_init(part);
	partition_z(part, n);
	mpf_set_z(res, part);
#if 0
static int last=100;
if (last < n) {
printf("duuude parti n=%d bits=%lu\n", n, mpz_sizeinbase(part, 2));
last = n;
}
#endif

	mpz_clear(part);
}

static double partition_big(double re_q, double im_q, int itermax, double param)
{
	// For radius to 600, we need maybe 50 decimals (340 bits)
	// For radius 6000, partition fun needs 120 decimals (400 bits)
	int nprec = 185;
	mpf_set_default_prec(nprec * 3.321);

	cpx_t sum, z; cpx_init(sum); cpx_init(z);
	mpf_t val; mpf_init(val);

	cpx_set_d(z, re_q, im_q);

	// cpx_exponential_genfunc(sum, z, nprec, partition);
	// cpx_exponential_genfunc(sum, z, nprec, foop);
	// cpx_exponential_genfunc_mpf(sum, z, nprec, parti_mpf);
	cpx_exponential_genfunc_mpf(sum, z, nprec, parti_z_mpf);
#define MAG 1
#if MAG

#if 0
	mpf_t gabs; mpf_init(gabs);
	cpx_abs(gabs, z);
	mpf_sqrt(gabs, gabs);
	mpf_neg(gabs, gabs);
	fp_exp(gabs, gabs, nprec);
	cpx_times_mpf(sum, sum, gabs);
	mpf_clear(gabs);
#endif

	cpx_abs(val, sum);

	double rv = mpf_get_d(val);

	// rv = log(1.0 + rv);
#endif

#if PHASE
	double rv = 0.5 + 0.5 * atan2(cpx_get_im(sum), cpx_get_re(sum))/M_PI;
#endif

	cpx_clear(sum);
	cpx_clear(z);
	mpf_clear(val);
	return rv;
}

// static double thue_morse_recip(long n) { return 1.0 / (1.0 + thue_morse(n)); }
static double thue_morse_rev(long n) { return 1.0 - thue_morse(n); }
static double xperiment(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	double complex z = re_q + I * im_q;
	double complex g = exp_genfunc_d(z, thue_morse_rev);
	return cabs(g);
}

static double necklace_mag(double re_q, double im_q, int itermax, double param)
{
	max_iter = itermax;
	long double complex z = re_q + I * im_q;
	// double complex g = exponential_genfunc(z, necklace);
	double complex g = ordinary_genfunc(z, necklace);
	return cabs(g);
	// return 0.5 + 0.5 * atan2(cimag(g), creal(g))/M_PI;
}

// ========================================================


__attribute__((constructor)) void decl_things() {
	DECL_HEIGHT("totient_ord_phase", totient_ord_phase);
	DECL_HEIGHT("totient_exp_phase", totient_exp_phase);
	DECL_HEIGHT("totient_exp_mag", totient_exp_mag);
	DECL_HEIGHT("totient_big", totient_big);
	DECL_HEIGHT("carmichael", carmichael);
	DECL_HEIGHT("carmichael_big", carmichael_big);
	DECL_HEIGHT("mobius_exp_mag", mobius_exp_mag);
	DECL_HEIGHT("mobius_big", mobius_big);
	DECL_HEIGHT("divisor_exp_mag", divisor_exp_mag);
	DECL_HEIGHT("divisor_big", divisor_big);
	DECL_HEIGHT("divisor_twist", divisor_twist);
	DECL_HEIGHT("divisor_uncircle", divisor_uncircle);
	DECL_HEIGHT("sigma_one", sigma_one);
	DECL_HEIGHT("sigma_two", sigma_two);
	DECL_HEIGHT("little_omega", little_omega_big);
	DECL_HEIGHT("unitary", unitary_big);
	DECL_HEIGHT("big_omega_exp_mag", big_omega_exp_mag);
	DECL_HEIGHT("big_omega_big", big_omega_big);
	DECL_HEIGHT("liouv_lambda", liouv_lambda);
	DECL_HEIGHT("mertens_m", mertens_m_exp_mag);

	DECL_HEIGHT("mangoldt_lambda", mangoldt_lambda_exp_mag);
	DECL_HEIGHT("mangoldt_lambda_big", mangoldt_lambda_big);
	DECL_HEIGHT("exp_mangoldt_lambda", exp_mangoldt_lambda_exp_mag);
	DECL_HEIGHT("exp_mango_big", exp_mango_big);
	DECL_HEIGHT("thue_morse", thue_morse_exp_mag);
	DECL_HEIGHT("thue_morse_big", thue_morse_big);
	DECL_HEIGHT("isqrt_big", isqrt_big);
	DECL_HEIGHT("partition_big", partition_big);
	DECL_HEIGHT("xperiment", xperiment);
	DECL_HEIGHT("necklace", necklace_mag);
}

// DECL_MAKE_HEIGHT(plot_big);
MAKE_HEIGHT;
