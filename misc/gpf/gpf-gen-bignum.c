/*
 * Generating functions for greatest prime factors.
 * Implementation in bignums.
 *
 * April 2016
 */

#include <math.h>
#include <stdio.h>

#include <mp-complex.h>

#include <gpf.h>
#include "gpf-gen-bignum.h"

/*
 * Ordinary generating function for the greatest common factor.
 */
void cpx_gpf_ordinary(cpx_t sum, cpx_t z, int prec)
{
	mpf_t zabs, gabs, epsi;
	mpf_init (gabs);
	mpf_init (zabs);
	mpf_init (epsi);
	mpf_set_ui(epsi, 1);
	mpf_div_2exp(epsi, epsi, (int)(3.321*prec));

	cpx_set_ui(sum, 0, 0);

	// falls apart if z is zero.
	cpx_abs(gabs, z);
	if (0 > mpf_cmp(gabs, epsi)) return;

	// Not defined for |z| > 1
	mpf_sub_ui(gabs, gabs, 1);
	if (0 > mpf_cmp(gabs, epsi)) return;

	cpx_t zn, term;
	cpx_init(zn);
	cpx_init(term);
	cpx_set(zn, z);

	for (int n=1; ; n++)
	{
		cpx_mul_ui(term, zn, gpf(n), 0);
		cpx_add(sum, sum, term);
		cpx_mul(zn, zn, z);

		cpx_abs(gabs, zn);
		mpf_mul_ui(gabs, gabs, n);

		cpx_abs(zabs, sum);
		mpf_div(gabs, gabs, zabs);

		// if (n * zn < epsi * sum) return;
		if (0 > mpf_cmp(gabs, epsi)) return;
	}
}

#if 0
/*
 * Exponential generating function for the greatest common factor.
 */
cpx_t cpx_gpf_exponential(long cpx_t x)
{
	long cpx_t sum = 0;
	long cpx_t xn = x;
	long double fact = 1.0;

	if (cabsl(x) < MAX_PREC) return x;

	for (int n=1; ; n++)
	{
		sum += gpf(n) * (xn * fact);
		xn *= x;
		fact /= n;
		if (n*cabsl(xn*fact) < MAX_PREC*cabsl(sum)) break;
		if (max_iter < n) break;
	}
scale = expl(2.0L * cabsl(x));
scale /= cabsl(x);
sum *= scale;
// printf("duuude %g sum=%g\n", cabs(x), cabs(sum));

	return sum;
}
#endif
