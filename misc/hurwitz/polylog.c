
/* polylog.c
 *
 * Implement Borwein-style polylogarithm 
 * on unit circle.
 *
 * Linas November 2006
 */

#include <math.h>
#include <stdio.h>
#include "binomial.h"

typedef struct {
	double re;
	double im;
} cplex;

static inline cplex cplex_zero (void)
{
	cplex z; z.re=0.0; z.im=0.0; return z;
}
static inline cplex cplex_one (void)
{
	cplex z; z.re=1.0; z.im=0.0; return z;
}

static inline cplex cplex_minus (cplex z)
{
	cplex rv;
	rv.re = -z.re;
	rv.im = -z.im;
	return rv;
}

static inline cplex cplex_add (cplex a, cplex b)
{
	cplex rv;
	rv.re = a.re + b.re;
	rv.im = a.im + b.im;
	return rv;
}

static inline cplex cplex_sub (cplex a, cplex b)
{
	cplex rv;
	rv.re = a.re - b.re;
	rv.im = a.im - b.im;
	return rv;
}

static inline cplex cplex_scale (double x, cplex z)
{
	cplex rv;
	rv.re = x * z.re;
	rv.im = x * z.im;
	return rv;
}

static inline cplex cplex_mult (cplex a, cplex b)
{
	cplex rv;
	rv.re = a.re * b.re - a.im * b.im;
	rv.im = a.re * b.im + b.re * a.im;
	return rv;
}

static inline cplex cplex_recip (cplex z)
{
	cplex rv;
	double mag = 1.0/ (z.re*z.re + z.im*z.im);
	rv.re = z.re * mag;
	rv.im = -z.im * mag;
	return rv;
}

static inline double cplex_modulus (cplex z)
{
	return sqrt (z.re*z.re + z.im*z.im);
}

static inline double cplex_phase (cplex z)
{
	return atan2 (z.im, z.re);
}

static inline cplex cplex_pow (cplex z, int n)
{
	cplex rv;
	if (50>n)
	{
		int k;
		rv = z;
		for (k=1; k<n; k++)
		{
			rv = cplex_mult (rv, z);
		}
	}
	else
	{
		fprintf (stderr, "pow=%d unimplemented\n", n);
	}
	return rv;
}

static inline cplex cplex_exp (cplex z)
{
	cplex rv;
	rv.re = rv.im = exp (z.re);
	rv.re *= cos (z.im);
	rv.im *= sin (z.im);
	return rv;
}

/* 
 * Return value of sum^k_{j=0} (n j) oz^j
 */
static cplex bee_k (int n, int k, cplex oz)
{
	int j;
	cplex pz = cplex_one();
	cplex acc = cplex_zero();

	double sgn = 1.0;
	for (j=0; j<=k; j++)
	{
		cplex term = pz;
		double bin = sgn * binomial (n,j);

		term = cplex_scale (bin, term);
		acc = cplex_add (acc, term);
		pz = cplex_mult(pz, oz);
		sgn = -sgn;
	}

	return acc;
}

/*
 * return estimate of polylog, Li_s(z) for estimator n.
 */
static cplex polylog_est (cplex s, cplex z, int n)
{
	int k;
	cplex oz = cplex_recip(z);

	cplex ska = z;
	ska.re -= 1.0;
	ska = cplex_mult (ska, cplex_recip(cplex_mult(z,z)));
	ska = cplex_pow (ska,n);
	
	cplex pz = cplex_one();
	cplex acc = cplex_zero();

	for (k=0; k<=2*n; k++)
	{
		/* the coefficient c_k of the polynomial */
		cplex ck = bee_k (n,k,oz);
		ck = cplex_sub (ck, ska);
		ck = cplex_mult (ck, pz);
		
		/* the inverse integer power */
		cplex dir = cplex_scale (log(k+1), s);
		dir = cplex_exp (cplex_minus(dir));

		/* put it together */
		cplex term = cplex_mult (ck, dir);
		acc = cplex_add (acc, term);

		pz = cplex_mult(pz, oz);
	}

	ska = cplex_recip(ska);
	acc = cplex_mult (acc, ska);
	acc = cplex_mult (acc, z);
	acc = cplex_minus (acc);

	return acc;
}

/* incomplete implementation */
static cplex hurwitz_beta (cplex s, double q)
{
	int n=41;
	cplex z;

	if (0.16 > q)
	{
		/* use the duplication formula to get into convergent region */
		cplex sm = cplex_minus(s);
		sm.re += 1.0;
		cplex ts = cplex_exp(cplex_scale (log(2), sm));
		cplex bt = hurwitz_beta (s, 2.0*q);
		bt = cplex_mult (bt, ts);

		cplex z = hurwitz_beta (s, q+0.5);
		z = cplex_sub(bt, z);
		return z;

	}
	else if (1.0e-10 < q)
	{
		/* normal case; within the convergence region */
		double ph = 2.0*M_PI*q;
		z.re = cos (ph);
		z.im = sin (ph);
		return polylog_est (s, z, n);
	}
	else 
	{
		return cplex_zero();
	}
}

main ()
{
	cplex z;

	cplex s;
	s.re = 2.0;
	s.im = 0.0;

	double q;
	for (q=0.0; q<1.0; q+=0.05)
	{
		z = hurwitz_beta (s, q);
		printf ("%g	%g	%g\n", q, z.re, z.im);
	}
}
