
/*
 * asymptote.c
 * Find saddle point of asympttic expansion
 *
 * Linas Vepstas January 2010
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>

int main()
{
	double complex a,b,c,d,h, xp, xn;
	double k;

	h = M_LN2 - I * M_PI;

	for (k=0.02; k<5.0; k+= 0.02)
	{
   	a = 1.0;
   	b = (k-1.0) * (1.0 - h) / h;
   	c = k * (2.0 - h) / h;

		d = b*b - 4.0 * a * c;
		d = csqrt(d);

		xp = 0.5 * (-b + d);
		xn = 0.5 * (-b - d);
		// printf("its k=%f %f+i%f %f+i%f\n", k, 
		//	creal(xp), cimag(xp), 
		//	creal(xn), cimag(xn));
		printf("%f	%f	%f\n", k, creal(xp), cimag(xp));
	}
}
