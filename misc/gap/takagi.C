
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

long double triangle (long double x)
{
	long double t = x - floorl(x);
	if (0.5L > t) return 2.0L*t;
	return 2.0L*(1.0L-t);
}

double square (double x)
{
	x -= floor (x);
	if (1.0 > x*2.0) return 0.0;
	return 1.0;
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

double takagi (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 1.0;
	for (k=0; k<1000000; k++)
	{
		acc += tw* triangle (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > tw) break;
	}

	return acc;
}

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

double dtakagi (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 4.0;
	for (k=2; k<10000000; k++)
	{
		double term = k* (k-1)*tw;
		acc += term* triangle (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > term) break;
	}

	return acc;
}

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

double takagi_prime (double w, double x)
{
	int k;
	double acc = 0.0;
	double tw = 1.0;
	double tp = 1.0;
	for (k=0; k<1000; k++)
	{
		acc += tw* square (tp*x);
		tp *= 2.0;
		tw *= w;
		if (1.0e-14 > tw) break;
	}

	return acc;
}

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

// ------------------------------------------------------------------
// ------------------------------------------------------------------

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

main (int argc, char *argv[])
{
	int i;

	// int nmax = 512;
	// int nmax = 431;
	int nmax = 3;

	if (argc <2)
	{
		printf ("Usage: %s <w-value>\n", argv[0]);
		exit (1);
	}
	double w = atof(argv[1]);

	for (i=1; i<nmax; i++)
	{
		double ts = 1.0;
		double x = i/((double)nmax);
		// double tw = takagi (w, x);
		// double ts = sin_takagi (w, x);
		// double tw = dtakagi (w, x);
		// double tw = log (takagi(w,x));
		// double tw = ttakagi(w,x);
		// double tw = div_takagi (w, x);
		// double tw = takagi_bumps (w, x);
		// double tw = takagi_prime (w, x);
		// double tw = takagi_exp (w, x);

		double tw = lytic (w, x);


		printf ("%d	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, tw, ts, tw-ts);
	}
}
