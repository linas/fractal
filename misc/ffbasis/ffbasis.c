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
 * A vector in the kernel of E aka a_k
 */
double kern(int k)
{
	int n;
	double sum = 0.0;

	double fourpi = - 4.0 * M_PI * M_PI;
	double numer = - 2.0 * M_PI;
	double fact = 1.0;
	for(n=0; n<130; n++)
	{
		double term = numer;
		double bino = binomial (2*n+k, 2*n);
		bino /= fact;

		term *= bino;
		sum += term;
		// printf("n=%d numer=%g  bino=%g term=%g sum=%g\n", n, numer, bino, term, sum);

		if (fabs(term) < 1.0e-30) break;

		numer *= fourpi;
		fact *= (2.0 * n + 3.0)  * (2.0*n + 2.0);
	}
	// printf(" ---------------------- above was k=%d\n", k);

	if (k%2 == 1) sum = - sum;
	return sum;
}

/**
 * The inverse of E
 */
double Einv_km(int k, int m)
{
	int j;
	if (1 == m && 0 == k) return 1.0;
	if (0 == k) return 0.0;

	double sum = 0.0;

	double pj = m*m;
	for (j=1; j<k; j++)
	{
		double term = kern(j);
		sum += pj * term;
		pj *= m;
	}

	sum *= (1.0+m) / pj;
	sum += kern(k);
	sum *= (1.0+m) / ((double) m);

	return sum;
}

void chk_Einv(void)
{
	int m,n, k;

	m=2;
	n=2;

	double sum = 0;
	for (k=0; k<130; k++)
	{
		double term = E_mk(m,k) * Einv_km(k,n);
		sum += term;
		printf("duude k=%d sum=%g term=%g\n", k, sum, term);
		if (k>10 && fabs(term) < 1.0e-10) break;
	}
}

/**
 * Check that kern is actually in the kernel of E
 * .. and it is, although rounding errors get really nasty around m=8
 */
void check_kern(int m)
{
	int k;

	printf("============ m= %d ==============\n", m);
	double sum = 0.0;
	for (k=0; k<1550; k++)
	{
		double term = E_mk(m, k) * kern(k);
		sum += term;
		printf("k=%d term=%g sum=%g\n", k, term, sum);

		if (fabs(term) < 1.0e-10 && k>10) break;
	}
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

void chk_log_everywhere()
{
	double x;

	for (x= 0.1; x < 15.2; x +=0.1123) {
		chk_log(x);
	}
}


void print_kern(void)
{
	int k;

	printf("#\n# The kernel series a_k\n#\n");
	for (k=0; k<248; k++)
	{
		double a_k = kern(k);
		printf("%d	%g\n", k, a_k);
	}
}

void check_kern_everywhere(void)
{
	int m;
	for (m=1; m<10; m++)
	{
		check_kern(m);
	}
}

int
main (int argc, char * argv[])
{
	//print_kern();

	chk_Einv();

	return 0;
}

