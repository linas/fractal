
/* polylog.c
 *
 * Implement Borwein-style polylogarithm 
 * on unit circle.
 *
 * Linas November 2006
 */

#include <math.h>
#include "binomial.h"

typedef struct {
	double re;
	double im;
} cplex;

static inline cplex cplex_one (void)
{
	cplex z; z.re=0.0; z.im=0.0; return z;
}
static inline cplex cplex_one (void)
{
	cplex z; z.re=1.0; z.im=0.0; return z;
}

static inline cplex cplex_invert (cplex z)
{
	cplex rv;
	double mag = 1.0/ (z.re*z.re + z.im*z.im):
	rv.re = z.re * mag;
	rv.im = -z.im * mag;
	return rv;
}

static inline cplex cplex_add (cplex a, cplex b)
{
	cplex rv;
	rv.re = a.re + b.re;
	rv.im = a.im * b.im;
	return rv;
}

static inline cplex cplex_scale (double x, cplex z);
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

/* 
 * Return value of sum^k_{j=0} (n j) oz^j
 */
cplex bee_k (int n, int k, cplex oz)
{
	int j;
	cplex pzre = cplex_one();
	cplex acc = cplex_zero();

	double sgn = 1.0;
	for (j=0; j<=k; j++)
	{
		cplex term = pxre;
		double bin = sgn * binomial (n,j);

		term = cplex_scale (bin, term);
		acc = cplex_add (acc, term);
		
		pzre = cplex_mult(pzre, oz);
		sgn = -sgn;
	}

	return acc;
}

main ()
{
	cplex z;
	z.re = -1.0;
	z.im = 0.0;

	z = bee_k (6,4, z);

	printf ("duude %g %g \n", z.re, z.im);
}
