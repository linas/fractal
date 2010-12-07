/*
 * newton-series.c
 *
 * Newton-series summation of the continued fraction
 * riemann integral. For S_12 only.
 *
 * Idea looks good on paper, but numeric check shows that the k and a1, a2 
 * cannot be interchanged.  Doing so makes the a1, a2 sums logarithmically 
 * divergent.
 *
 * Linas Vepstas December 2010
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

complex alpha(complex ess, complex *betarr, complex *gamarr,  int kmax, int a1max, int a2max) 
{
	unsigned int k, na1, na2;
	complex alpha_sum = 0.0;

	for (na1=1; na1<a1max; na1++)
	{
		for (na2=1; na2<a2max; na2++)
		{
			if (na1 == na2) continue;  // c==0, d==0

			double a1 = na1;
			double a2 = na2;
			double b = a1 - a2;
			double c = - (a1 * a2 + 1.0) * b;
			double d = 1.0 + a2 * b;
			double xlo = a2 / (1.0 + a1 * a2);
			double xhi = (1.0 + a2) / (1.0 + a1 + a1 * a2);
			double greb = - d / c;
			complex scale = cpow(greb, ess-1.0) / (c*c);

			complex term = log((xhi+greb) /(xlo+greb));
			term *= scale;
			alpha_sum += term;

			double ced = c/d;
			double exblo = 1.0 + ced*xlo;
			double exbhi = 1.0 + ced*xhi;

			double klo = exblo;
			double khi = exbhi;
			for (k=1; k<kmax; k++)
			{
				term = khi - klo;
				term *= scale;

				betarr[k] += term;
				gamarr[k] += term*b*c;

				klo *= exblo;
				khi *= exbhi;
			}
		}
	}
	return alpha_sum;
}

int main (int argc, char * argv[])
{
	int k;
	int amax = atoi(argv[1]);
	k = atoi(argv[2]);

	complex ess = 0.5 + I*7.0;

	int kmax = k+1;
	complex *betarr = (complex *) malloc(kmax * sizeof(complex)); 
	complex *gamarr = (complex *) malloc(kmax * sizeof(complex)); 

	double last = 0;
	while (1) {
	complex yo = alpha(ess, betarr, gamarr, kmax, amax, amax);

	double re = creal (yo);
	double im = cimag(yo);

	// printf("duude %d %20.18g %20.18g\n", amax, re, im);
	//for (k=1; k<kmax; k++)
	{
		re = creal (betarr[k]);
		im = cimag (betarr[k]);

		printf("%d	k=%d %20.18g %20.18g delta=%g\n", amax, k, re, im, re-last);
		last = re;
	}
	amax *= 2;
	}
	return 0;
}
