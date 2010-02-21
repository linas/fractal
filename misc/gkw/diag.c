
/*
 * diag.c
 * check saddle point
 *
 * Linas Vepstas January 2010
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>

int main()
{
	double complex h, x02, x0, t, a;

	h = M_LN2 - I * M_PI;
	x02 = 1 - 2.0/h;
	x0 = csqrt(x02);

	printf("x_0 = %f +i %f\n", creal(x0), cimag(x0));

	t = (1.0-x0) / (1.0+x0);
	printf("1-x/1+x = %f +i %f\n", creal(t), cimag(t));

	a = h * x0 + clog(t);
	printf("a = %f +i %f\n", creal(a), cimag(a));
}
