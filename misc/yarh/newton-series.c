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
 * The "obvious" work-around to this is to regulate the sums.
 *
 * Linas Vepstas December 2010
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

complex alpha(complex ess, double tee, complex *betarr, complex *gamarr,  int kmax) 
{
	unsigned int k, na1, na2;
	complex alpha_sum = 0.0;

	double reg = exp(-tee);
	double rega1 = 1.0;

#define EPSILON 1.0e-10
	for (na1=1; 1; na1++)
	{
		rega1 *= reg;
		// Terminate according to the regulator.
		if (rega1 < EPSILON) goto done;

		double rega2 = 1.0;
		for (na2=1; 1; na2++)
		{
			rega2 *= reg;
			if (na1 == na2) continue;  // c==0, d==0

			// Terminate according to the regulator.
			if (rega1*rega2 < EPSILON) break;

			double a1 = na1;
			double a2 = na2;
			double b = a1 - a2;
			double c = - (a1 * a2 + 1.0) * b;
			double d = 1.0 + a2 * b;
			double xlo = a2 / (1.0 + a1 * a2);
			double xhi = (1.0 + a2) / (1.0 + a1 + a1 * a2);
			double greb = - d / c;
			complex scale = rega1*rega2 * cpow(greb, ess-1.0) / (c*c);

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

done:
	return alpha_sum;
}

int main (int argc, char * argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s <kmax> <tee>\n", argv[0]);
		exit(1);
	}

	int kmax;
	kmax = atoi(argv[1]);
	double tee = atof(argv[2]);

	complex ess = 0.5 + I*7.0;

	complex *betarr = (complex *) malloc(kmax * sizeof(complex)); 
	complex *gamarr = (complex *) malloc(kmax * sizeof(complex)); 

	complex yo = alpha(ess, tee, betarr, gamarr, kmax);

	double re = creal (yo);
	double im = cimag(yo);
	printf("tee=%g %20.18g %20.18g\n", tee, re, im);

	int k;
	for (k=1; k<kmax; k++)
	{
		re = creal (betarr[k]);
		im = cimag (betarr[k]);

		printf("k=%d %20.18g %20.18g\n", k, re, im);
	}
	return 0;
}
