/*
 * euler.c
 *
 * Implement Euler summation for some complex-valued arithmetic function.
 * viz, given f(n) for 0 < n return
 *
 * sum_{n=0}^\infty 2^{-(n+1)} \sum_{k=0}^n {n choose k} f(k+1)
 *
 * See wikipedia for more about Euler summation.
 *
 * March 2019  -- Linas Vepstas
 */

#include <stdio.h>
#include "binomial.h"
#include "euler.h"

complex euler_sum(arithmetic fun)
{
	complex sum = 0.0;
	long double tn = 0.5;

	// DBL_MAX = 1.79769e+308
	// so log_2 of that is approx 1024
	int n = 0;
	for (; n<1000; n++)
	{
		complex term = 0.0;
		for (int k=0; k<=n; k++)
		{
			term += binomial(n, k) * fun(k+1);
		}
		term *= tn;
		sum += term;
if (cabs(term) < 1.0e-20) {printf("duuuude n=%d tn=%Lg\n", n, tn);}
		if (cabs(term) < 1.0e-20) return sum;
		tn *= 0.5;
	}
	fprintf(stderr, "Warning: Euler-sum double overflow!\n");

	// Keep going some more. Watch out for exponent overflow.
	// LDBL_MAX = 1.18973e+4932
	// so log_2 of that is approx 16384
	for (; n<16300; n++)
	{
		complex term = 0.0;
		for (int k=0; k<=n; k++)
		{
			// Avoid insane exponents
			long double bino = binomial(n, k) * tn;
			term += bino * fun(k+1);
		}
		sum += term;
printf("its %d sum=%f+i%f term=%g+i%g tn=%Lg\n", n, creal(sum), cimag(sum),
creal(term), cimag(term), tn);
if (cabs(term) < 1.0e-20) {printf("duuuude n=%d tn=%Lg\n", n, tn);}
		if (cabs(term) < 1.0e-20) return sum;
		tn *= 0.5;
	}
	fprintf(stderr, "Error: Euler-sum non-convergence!\n");
	return sum;
}

#define TEST
#ifdef TEST

#include <float.h>
#include <stdio.h>

complex z = 0.0;
complex zee(int n) { return cpow(z, n-1); }

int main()
{
	printf("DBL_MAX=%g LDBL_MAX =%Lg\n", DBL_MAX, LDBL_MAX);
	// print_binomial();
	z = 0.7 + I*0.5;
	z = 0.975;
	complex rslt = euler_sum(zee);
	complex xpct = 1.0 / (1.0-z);
	printf("z=%f + i%f  got %f + i%f expect %f + i%f\n",
		creal(z), cimag(z), creal(rslt), cimag(rslt), creal(xpct), cimag(xpct));
	complex diff = rslt - xpct;
	printf("diff %g + i%g\n", creal(diff), cimag(diff));
}
#endif
