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
#ifndef EULER_H__
#define EULER_H__

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

//  C and C++ is fugnuts insane in complex support.
#define complex _Complex

typedef complex arithmetic(int);

complex euler_sum(arithmetic);


#ifdef __cplusplus
};
#endif

#endif /* EULER_H__ */

#include "binomial.h"

complex euler_sum(arithmetic fun)
{
	complex sum = 0.0;
	double tn = 0.5;
	int n = 0;
	while (1)
	{
		complex term = 0.0;
		for (int k=0; k<=n; k++)
		{
			term += binomial(n, k) * fun(k+1);
		}
		term *= tn;
		sum += term;
		if (cabs(term) < 1.0e-20) return sum;
		n++;
		tn *= 0.5;
	}
	return 0.0;
}

#define TEST
#ifdef TEST

#include <stdio.h>

complex z = 0.0;
complex zee(int n) { return cpow(z, n-1); }

int main()
{
	z = 0.7 + I*0.5;
	complex rslt = euler_sum(zee);
	complex xpct = 1.0 / (1.0-z);
	printf("z=%f + i%f  got %f + i%f expect %f + i%f\n",
		creal(z), cimag(z), creal(rslt), cimag(rslt), creal(xpct), cimag(xpct));
	complex diff = rslt - xpct;
	printf("diff %g + i%g\n", creal(diff), cimag(diff));
}
#endif
