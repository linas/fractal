
/* 
 * takagi.C
 *
 * draw the Takagi curve
 *
 * wow! checkout tagaki bumps!  wonder if one could make a sin(pi/x) 
 * approximation out of this...
 *
 * Linas October 2004
 */

#include <complex.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_sf_zeta.h>
#include "Farey.h"

ContinuedFraction far;

long double triangle (long double x)
{
	long double t = x - floorl(x);
	if (0.5L > t) return 2.0L*t;
	return 2.0L*(1.0L-t);
}

double square (double x)
{
	x -= floor (x);
	if (1.0 > x*2.0) return -2.0;
	return 2.0;
}

long double parabola_down (long double x)
{
	long double t = x - floorl(x);
	return 4.0L*t*(1.0L-t);
}

long double parabola_up (long double x)
{
	long double t = x - floorl(x);
	if (0.5 < t) t = 1.0-t;
	return 4.0L*t*t;
}

double bumps (double x)
{
	x -= floor (x);
	if (1.0 > x*3.0) return 0.0;
	if (2.0 < x*3.0) return 0.0;
	return triangle (3.0*x);
}

double saw (double x)
{
	x *= 3.0;
	if (1.0 > x) return 0.0;
	if (2.0 < x) return 0.0;
	return 2.0-x;
}

double mink(double x)
{
	x = triangle (x);
	far.SetReal (x);
	x = far.ToFarey();
	return x;
}

long double sq (long double x)
{
	long double t = x - floorl(x);
	if (1.0L > 3.0L*t) return 2.0*t/(1.0-t);
	return 0.5*(1.0-t)/t;;
}

// the farey/isola map
long double pointy (long double x)
{
	long double t = x - floorl(x);
	if (0.5L < t) return (1.0L-t)/t;
	return t/(1.0L-t);
}


/* The main, core basic takagi curve */
long double takagi (long double w, long double x)
{
	int k;
	long double acc = 0.0L;
	long double tw = 1.0L;
	long double tp = 1.0L;
	for (k=0; k<50; k++)
	{
		long double term = tw* triangle (tp*x);
		// long double term = tw* parabola_down (tp*x);
		// long double term = tw* parabola_up (tp*x);
		acc += term;
		tp *= 2.0L;
		tw *= w;
		if (1.0e-16 > tw) break;
	}

	return acc;
}

/* The main, core basic isola curve */
long double isola (long double w, long double x)
{
	int k;
	long double acc = 0.0L;
	long double tw = 1.0L;
	long double xit = x;
	for (k=0; k<50; k++)
	{
		xit = pointy (xit);
		long double term = tw * xit;
		acc += term;
		tw *= w;
		if (1.0e-16 > tw) break;
	}

	return acc;
}

long double iter_tak (long double w, long double x)
{
	int k;
	long double acc = 0.0L;
	long double tw = 1.0L;
	long double xit = x;
	for (k=0; k<50; k++)
	{
		// xit = triangle (xit);
		// xit = parabola_up (xit);
		xit = sq (xit);
		// xit = parabola_down (xit);
		long double term = tw * xit;
		acc += term;
		tw *= w;
		if (1.0e-16 > tw) break;
	}

	return acc;
}

/* The main, core basic takagi curve */
// XXX the triangle is accureate only to 50 bits or so,
// so any further and we need arbit-precision.
// which means that this func is not accurate below s=1.5 or so.
long double dirichlet_takagi (long double s, long double x)
{
	int k;
	long double acc = 0.0L;
	long double tp = 1.0L;
	for (k=1; k<50; k++)
	{
		long double tw = k;
		tw = powl (tw, -s);
		long double term = tw* triangle (tp*x);
		// long double term = tw* parabola_down (tp*x);
		// long double term = tw* parabola_up (tp*x);
		acc += term;
		tp *= 2.0L;
		if (1.0e-16 > tw) break;
	}
printf ("duude last term=%d  %g\n", k, pow (k, -s));

	return acc;
}

/* a stupid, uninteresting 2*l+1 variation *.
long double takagi_l (long double w, int l, long double x)
{
	int k;
	long double acc = 0.0L;
	long double tw = 1.0L;
	long double tp = 1.0L;
	for (k=0; k<100; k++)
	{
		long double term = tw* triangle ((2*l+1)*tp*x);
		acc += term;
		tp *= 2.0L;
		tw *= w;
		if (1.0e-16 > tw) break;
	}

	return acc;
}

/* crude conjagte attempt */
double sin_takagi (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 1.0;
	for (k=0; k<10000000; k++)
	{
		acc += tw* (1.0 - cos(tp*x*(2.0*M_PI)))*0.5;
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > tw) break;
	}

	return acc;
}

/* derivative of takagi with respect to w */
double dtakagi (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 2.0;
	for (k=1; k<1000; k++)
	{
		double term = k*tw;
		acc += term* triangle (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > term) break;
	}

	return acc;
}

/* stupid, dorky change-of-variable on derivative of takagi */
double ttakagi (double s, double x)
{
	int k;
	double exps = exp(-s);
	double w = exps;
printf ("# duude w=%g\n", w);
	double acc = 0.0;
	double tw = 1.0;
	double tp = 2.0;
	for (k=1; k<10000000; k++)
	{
		acc += k* tw* triangle (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > tw) break;
	}

	return acc*exps;
}

/* gauusian regulated takagi */
double div_takagi (double cut, double x)
{
	int k;
	double w = 2.0;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 1.0;
	for (k=0; k<100000000; k++)
	{
		double reg = exp (-k*k*cut*cut);
		double term = reg * tw;
		acc += term * triangle (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > term) break;
	}

	acc *= 2.0 *cut*cut;

	return acc;
}

double takagi_bumps (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 1.0;
	for (k=0; k<1000; k++)
	{
		acc += tw* bumps (tp*x);
		// acc += tw* saw (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > tw) break;
	}

	return acc;
}

/* derivative of takagi with respect to x */
double takagi_prime (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 1.0;
	for (k=0; k<1000; k++)
	{
		acc += tw* tp* square (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-16 > tw) break;
	}

	return acc;
}

/* entire-function variant of takagi, done by stiking in factorial */
double takagi_exp (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 1.0;
	double fact = 1.0;
	for (k=0; k<1000; k++)
	{
		double term = tw / fact;
		acc += term* triangle (tp*x);
		tp *= 2.0;
		tw *= w;
		fact *= k+1;
		if (1.0e-14 > term) break;
	}

	return acc;
}

double takagi_mink (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 1.0;
	for (k=0; k<32; k++)
	{
		double term = tw* mink (tp*x);
		acc += term;
// printf ("# k=%d term=%g acc=%g\n", k, term, acc);
// fflush (stdout);
		tp *= 2.0;
		tw *= w;
		if (1.0e-16 > tw) break;
	}

	return acc;
}

#if FAILED_ANALYTIC_CONTINUATION
// ------------------------------------------------------------------
// ------------------------------------------------------------------
// failed attempt at brute-force analytic continuation

// ======================================================
// brute-force factorial function
inline long double factorial (int n)
{
	int k;
	long double fac = 1.0L;
	for (k=2; k<=n; k++)
	{
		fac *= (long double) k;
	}
	if (0>n) fac = 0.0L;
	return fac;
}

// ======================================================
// brute-force binomial coefficent
// must have m<=n, m>=0
inline long double binomial (int n, int m)
{
	if (0>m) return 0.0L;
	if (2*m < n) m = n-m;
	int l = n-m;
	if (0>l) return 0.0L;
	int k;
	long double bin = 1.0L;
	for (k=1; k<=l; k++)
	{
		bin *= (long double) (m+k);
	}
	bin /= factorial (l);
	return bin;
}

/* These complex functions do analytic continuation of takagi to
 * w greater than one (actually between 1 and 3) */
#define Complex complex<long double>

#define EPS 1.0e-8

Complex trin (int n, long double x)
{
	int k;
	long double twon = 1.0L;
	if (48 < n) return 0.0L;

	for (k=0; k<n; k++)
	{
		twon *= 2.0L;
	}
	Complex term = triangle (twon * x);
	// printf ("duude trianlge (%d, %g) = %g\n", n, x, real(term));
	return term;
}

Complex do_coeff(Complex a, int k, Complex (*fn)(int, long double), long double x)
{
	int n;
	Complex acc = 0.0L;
	Complex an = 1.0L;
	long double cut = 1.0/20.0;  // XXX need to extrapolate to zero 
	for (n=0; n<10000; n++)
	{
		Complex term = binomial(n+k,k);
		long double reg = expl (-cut*cut*n*n);
		term *= reg;
		term *= an * fn(n+k,x);
		acc += term;
		an *= a;

// printf ("n=%d term=%Lg + i %Lg   acc = %Lg + i %Lg\n", n, real(term), imag(term), real (acc), imag (acc));
		if (EPS > abs(term)) break;
	}
// printf ("------\n\n");
	return acc;
}

Complex a_coeff(int k, long double x)
{
	Complex a;
	a = Complex(0.5L, 0.5L);

// printf ("start acoeff k=%d x=%Lg\n", k, x);
	Complex term = do_coeff (a, k, trin, x);
printf ("end acoeff k=%d x=%Lg coeff=%Lg + i %Lg\n", k, x, real(term), imag (term));
	return term;
}

Complex b_coeff(int k, long double x)
{
	Complex a;
	a = 0.5L;

printf ("start bcoeff k=%d x=%Lg\n", k, x);
	Complex term = do_coeff (a, k, a_coeff, x);
printf ("end bcoeff k=%d x=%Lg coeff=%Lg + i %Lg\n", k, x, real(term), imag (term));
	return term;
}

Complex c_coeff(int k, long double x)
{
	Complex a;
	a = 0.5L;

printf ("start ccoeff k=%d x=%Lg\n", k, x);
	Complex term = do_coeff (a, k, b_coeff, x);
printf ("end ccoeff k=%d x=%Lg coeff=%Lg + i %Lg\n", k, x, real(term), imag (term));
	return term;
}

Complex d_coeff(int k, long double x)
{
	Complex a;
	a = Complex(0.5L, -0.5L);

printf ("start dcoeff k=%d x=%Lg\n", k, x);
	Complex term = do_coeff (a, k, c_coeff, x);
printf ("end dcoeff k=%d x=%Lg coeff=%Lg + i %Lg\n", k, x, real(term), imag (term));
	return term;
}

long double lytic (long double w, long double x)
{
	int k;
	w -= 2.0L;

	long double acc = 0.0;
	long double wk = 1.0;
	for (k=0; k<10000; k++)
	{
		Complex cterm = d_coeff (k, x);
		double term = real (cterm);
		term *= wk;
		acc += term;
		if (EPS > fabs(term)) break;
		wk *= w;
	}

	return acc;
}
#endif
// ===============================================
// ===============================================

main (int argc, char *argv[])
{
	int i;

	// int nmax = 512;
	// int nmax = 2048;
	int nmax = 12;

	if (argc <2)
	{
		printf ("Usage: %s <w-value>\n", argv[0]);
		exit (1);
	}
	double w = atof(argv[1]);

	for (i=0; i<nmax; i++)
	{
		double x = i/((double)nmax);
		// double ts = isola (w, x);
		// double tw = takagi (w, x);
		// double tw = iter_tak (w, x);

		double tw = dirichlet_takagi (w, x);
		double ts = 2.0*gsl_sf_zeta (w) / 3.0 -0.5 - pow(2.0, -w)/3.0;

		// double ts = sin_takagi (w, x);
		// double tw = dtakagi (w, x);
		// double tw = log (takagi(w,x));
		// double tw = ttakagi(w,x);
		// double tw = div_takagi (w, x);
		// double tw = takagi_bumps (w, x);
		// double tw = takagi_prime (w, x);

		// double tw = lytic (w, x);

#ifdef HALF_SYM
		// g
		double tw = takagi (w, 0.5*x);
		double ts = x + w*takagi (w, x);

		// g*g
		double tw = takagi (w, 0.25*x);
		double ts = 0.5*x + w*x + w*w*takagi (w, x);

		// g*g*g
		double tw = takagi (w, 0.125*x);
		double ts =  x*(0.25+0.5*w +w*w) +w*w*w*takagi (w,x);

		// r*g*g*g
		double tw = takagi (w, 1.0-0.125*x);
		double ts = x*(0.25+0.5*w +w*w) +w*w*w*takagi (w,x);

		// g^4
		double tw = takagi (w, 0.0625*x);
		double ts = 0.125*x*(1.0-16.0*w*w*w*w)/(1.0-2.0*w) + w*w*w*w*takagi (w, x);
		// gr
		double tw = takagi (w, 0.5-0.5*x);
		double ts = 1.0 - x + w*takagi (w,x);

		// grg
		double tw = takagi (w, 0.5-0.25*x);
		double ts = 1.0 +x*(-0.5+w) +w*w*takagi (w, x);

		// grg^3
		double tw = takagi (w, 0.5-0.0625*x);
		double ts = 1.0 + x*(-0.125+0.25*w +0.5*w*w + w*w*w) +w*w*w*w*takagi (w,x);

		double tw = takagi (w, 0.5*(1.0+x));
		double ts = 1.0-x + w*takagi (w, x);
		double tw = takagi (w, 1.0+0.0625*(x-1.0));
		double ts = 0.125*(1.0-x)*(1.0-16.0*w*w*w*w)/(1.0-2.0*w) + w*w*w*w*takagi (w, x);
		double tw = takagi (w, 0.75-0.25*x);
		double ts = 0.5 -2.0*w*w + x*(0.5-2.0*w +2.0*w*w);
		ts /= 1.0-2.0*w;
		ts += w*w*takagi (w, x);

		double tw = takagi_l (w, 1, 0.5*x);
		double ts = takagi_l (w, 0, 3.0*x/2.0);
		double tw = takagi_l (w, 2, 0.5*x);
		double ts = takagi_l (w, 1, 5.0*x/6.0);
		double tw = takagi_l (w, 3, 0.5*x);
		double ts = takagi_l (w, 2, 7.0*x/10.0);
		int k=5;
		double tw = takagi_l (w, k, 0.5*x);
		double ts = takagi_l (w, k-1, (2*k+1)*x/(2*(2*k-1)));

		int k=4;
		double tw = takagi_l (w, k, 0.5*x);
		double ts = (2*k+1)*x + w*takagi_l (w, k, x);
#endif
		// double tw = takagi_mink (w, x);
		// double ts = 1.0;

		// double tw = takagi (w, x);
		// double ts = takagi_prime (0.5*w, x);

		// double tw = takagi (w, 0.5-0.0625*x);
		// double ts = 0.125 + 0.25*w + 0.5*w*w - x*(0.125+0.25*w* +0.5*w*w - w*w*w);
		// ts += w*w*w*w*takagi (w, x);

		// double tw = takagi (w, 0.125*x);
		// double ts =  x*(0.25+0.5*w +w*w) +w*w*w*takagi (w,x);

#ifdef FOUR_DIM_REP
		// g
		double tw = takagi (w, 0.5*x);
		double ts = 2.0*x - x*x + w*takagi (w, x);

		// g*g
		double tw = takagi (w, 0.25*x);
		double ts = x*(1.0+2.0*w) - x*x*(w+0.25) + w*w*takagi (w, x);

		// g*r
		double tw = takagi (w, 0.5*(1.0-x));
		double ts = 1.0 -x*x + w*takagi (w, x);

		// r * g
		double tw = takagi (w, (1.0-0.5*x));
		double ts = 2.0*x -x*x + w*takagi (w, x);

		// grg^2
		double tw = takagi (w, (0.5-0.125*x));
		double ts = 1.0 + x*(w+2.0*w*w) - x*x*(w*w+0.25*w+0.0625) + w*w*w*takagi (w, x);

#endif
		// double tw = triangle (128.0*(0.5*x));
		// tw += triangle (128.0*(0.5+0.5*x));
		// tw *= 0.5;
		// double ts = triangle (64.0*x);

		printf ("%d	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, tw, ts, tw-ts);
		fflush (stdout);
	}
}
