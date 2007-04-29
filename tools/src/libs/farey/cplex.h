
/** 
 * cplex.h
 *
 * Implement an actual, usable complex datatype,
 * unlike the ANSI C++ horseshit.
 *
 * Linas November 2006
 */

#include <math.h>
#include <stdio.h>

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
static inline cplex cplex_set (double x, double y)
{
	cplex z; z.re=x; z.im=y; return z;
}

static inline cplex cplex_conj (const cplex z)
{
	cplex zb; zb.re=z.re; zb.im=-z.im; return zb;
}

static inline cplex cplex_times_i (const cplex z)
{
	cplex rv;
	rv.re = -z.im;
	rv.im = z.re;
	return rv;
}

static inline cplex cplex_neg (const cplex z)
{
	cplex rv;
	rv.re = -z.re;
	rv.im = -z.im;
	return rv;
}

static inline double cplex_mag (const cplex z)
{
	return sqrt (z.re*z.re + z.im*z.im);
}

static inline cplex cplex_add (const cplex a, const cplex b)
{
	cplex rv;
	rv.re = a.re + b.re;
	rv.im = a.im + b.im;
	return rv;
}

static inline cplex cplex_sub (const cplex a, const cplex b)
{
	cplex rv;
	rv.re = a.re - b.re;
	rv.im = a.im - b.im;
	return rv;
}

static inline double cplex_dist (const cplex a, const cplex b)
{
	return cplex_mag (cplex_sub(a,b));
}

static inline cplex cplex_scale (double x, const cplex z)
{
	cplex rv;
	rv.re = x * z.re;
	rv.im = x * z.im;
	return rv;
}

static inline cplex cplex_mul (const cplex a, const cplex b)
{
	cplex rv;
	rv.re = a.re * b.re - a.im * b.im;
	rv.im = a.re * b.im + b.re * a.im;
	return rv;
}

static inline cplex cplex_recip (const cplex z)
{
	cplex rv;
	double mag = 1.0/ (z.re*z.re + z.im*z.im);
	rv.re = z.re * mag;
	rv.im = -z.im * mag;
	return rv;
}

static inline cplex cplex_div (const cplex a, const cplex b)
{
	cplex rv, deno;
	deno = cplex_recip (b);
	rv = cplex_mul (a, deno);
	return rv;
}

static inline double cplex_modulus (const cplex z)
{
	return sqrt (z.re*z.re + z.im*z.im);
}

static inline double cplex_phase (cplex z)
{
	return atan2 (z.im, z.re);
}

static inline cplex cplex_exp_itheta (double theta)
{
	cplex z;
	z.re = cos(theta);
	z.im = sin(theta);
	return z;
}

/** return z^n */
static inline cplex cplex_pow_ui (const cplex z, unsigned int n)
{
	cplex rv;
	rv = cplex_one ();
	if (350>n)
	{
		int k;
		for (k=0; k<n; k++)
		{
			rv = cplex_mul (rv, z);
		}
	}
	else
	{
		fprintf (stderr, "pow=%d unimplemented\n", n);
	}
	return rv;
}

static inline cplex cplex_exp (const cplex z)
{
	cplex rv;
	rv.re = rv.im = exp (z.re);
	rv.re *= cos (z.im);
	rv.im *= sin (z.im);
	return rv;
}

/** return n^s */
static inline cplex cplex_ui_pow (unsigned int n, cplex s)
{
	double lnn = log (n);
	s.re *= lnn;
	s.im *= lnn;
	return cplex_exp (s);
}

/** return q^s */
static inline cplex cplex_d_pow (double q, cplex s)
{
	double lnn = log (q);
	s.re *= lnn;
	s.im *= lnn;
	return cplex_exp (s);
}

