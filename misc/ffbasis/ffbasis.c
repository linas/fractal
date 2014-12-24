/**
 * ffbasis.c
 *
 * Fallig factorial basis for the divisor function
 *
 * Linas Vepstas Decmber 2014
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "binomial.h"

/**
 * The operator that transforms binomial (falling factorial) into a
 * Dirichlet series.
 *
 * We assume m>1 and that k >=0
 */
double E_mk(int m, int k)
{
	double em = m;
	em = (1.0 - em) / em;
	em = pow(em, ((double) k));
	return em;
} 

/**
 * The inverse of E
 */
dboule Einv_km(int k, int m)
{
}

/**
 * Consider the binomial b(s,k) = s! / (s-k)! k! for s complex
 * and k integer.  Verify that the relation 
 *   1/m^s = sum_k=0^\infty E_mk b(s,k)
 * holds.  (and indeed it soes seem to for all m and s.)
 * arguments are s and m.
 * Print a failure message if the relation fails to hold.
 */
void chk_E(long double complex s, int m)
{
	int k;
	long double complex sum = 0.0;
	for (k=0; k<800; k++)
	{
		double emk = E_mk(m, k);
		long double complex bin = cbinomial(s, k);
		sum += emk * bin;
	}

	long double complex psi = cpowl(m, -s);

	long double complex diff = sum - psi;

	// printf("its m=%d sum=%lf vs psi=%lf \n", m, creal(sum), creal(psi));
	// printf("its m=%d diff= %Lg + i %Lg  \n", m, creall(diff), cimagl(diff));
	long double adiff = cabsl(diff);

	if (1.0e-6 < adiff)
		printf("Fail at s=%f +i %f m=%d diff= %Lg \n", creal(s), cimag(s), m, adiff);
}

/**
 * call chk_E for a variety of points on the complex plane.
 */
void chk_E_everywhere(void)
{
	int m;
	double x,y;

#define MAX 12
	for (m=1; m<MAX ; m++)
	{
		for (x= - 5.0; x < 5.0; x +=0.30345) {
			for (y= - 5.0; y < 5.0; y += 0.4089345) {
				long double complex s = x + y *I;
				chk_E(s, m);
			}
				printf("."); fflush (stdout);
		}
		printf("\ndone with m=%d\n", m);
	} 
}

/**
 * Check the well-known formula
 *   log(x) = sum_k=1^\infty x^k / k
 * (It shows up in E-related calculations.)
 */
void chk_log(double x)
{
	int k;
	double y = (x-1.0) / x;
	double yp = y;
	double sum = 0.0;
	for (k=1; k<153; k++)
	{
		sum += yp / k;
		yp *= y;
	}

	printf("x=%f log(x)=%f  diff=%g\n", x, log(x), log(x)-sum);
}

int
main (int argc, char * argv[])
{
	double x;

	for (x= 0.1; x < 15.2; x +=0.1123) {
		chk_log(x);
	}
	return 0;
}

