/*
 * exact-series.c
 *
 * Solution via series
 *
 * Linas Vepstas December 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "binomial.h"

/* Integral of (x^(s-1))/(x+alpha) dx
 */
long double complex eff (long double complex s, long double alpha, long double x)
{
	int k;
	long double complex term, sum;
	long double opxa = 1.0L + x/alpha;
	long double opxak = opxa;

	sum = 0.0L;

	for (k=1; k<155456123; k++)
	{
		term = cbinomial(s,k);
		// printf("k=%d bino =%g +I %g\n", k, creal(term), cimag(term));
		term *= opxak;
		term /= (long double) k;
		if (k%2)
		{
			sum -= term;
		}
		else
		{
			sum += term;
		}

		// printf("k=%d term=%g +I %g\n", k, creal(term), cimag(term));
		double tm = cabs(term);
		if (tm < 1.0e-30) break;

		opxak *= opxa;
	}
	// printf("\nsum =%g +I %g\n", creal(sum), cimag(sum));

	long double xa = x+alpha;
	if (xa == 0.0L)
	{
		fprintf(stderr, "Error: unexpected sign for x+alpha=%Lg x=%Lg alpha=%Lg\n", xa, x, alpha);
		exit(1);
	}
	else if (xa > 0.0L)
	{
		term = logl(xa);
		sum += term;
	}
	else
	{
		term = logl(-xa);
		sum -= term;
	}
	// printf("sumpost log =%g +I %g\n", creal(sum), cimag(sum));

	/* Scale by (-alpha)^s */
	if (0.0L < alpha)
	{
		fprintf(stderr, "Oh no Mr. Bill, unexpected alpah=%Lg\n", alpha);
		exit(1);
	}
	term = s * logl(-alpha);
	term = cexpl(term);

	sum *= term;

	// printf("sumpost norm =%g +I %g\n\n", creal(sum), cimag(sum));
	return sum;
}

/* Integral of (x^s) dx = (x^(s+1))/(s+1)
 */
long double complex effo (long double complex s, long double x)
{
	long double complex term = logl(x);
	term *= (s+1.0L);
	term = cexpl(term);
	term /= (s+1.0L);

	return term;
}

/* Integral of (x^(s-1))/(x+alpha) dx
 * however, converges differently
 */
long double complex special_eff (long double complex s, long double alpha, long double x)
{
	int n;
	long double complex term, sum;
	long double rat, nrat;

	rat = alpha / x;
	nrat =1.0L;

	sum = 0.0L;
	for (n=0; n<57123456; n++)
	{
		term = nrat / (s-1.0L - (long double) n);
		sum += term;

		double tm = cabs(term);
		if (tm < 1.0e-30) break;

		nrat *= rat;
	}

	term = logl(x);
	term *= (s-1.0L);
	term = cexpl(term);
	sum *= term;

	return sum;
}

/* Integral of s_12 */

long double complex 
gral_s12(long double complex s, double epsi,
         unsigned int a1max, double *termerr)
{
	unsigned int na1, na2;
	long double complex term;
	long double complex sum = 0.0L;
	double cnt = 0.0;
	
	for (na1=1; na1<a1max; na1++)
	{
		for (na2=1; 1; na2++)
		{
			long double a1 = na1;
			long double a2 = na2;
			long double xlo = a2 / (1.0L + a1 * a2);
			long double xhi = (1.0L + a2) / (1.0L + a1 + a1 * a2);

			long double b = a1 - a2;
			long double c = - (a1 * a2 + 1.0L) * b;
			long double d = 1.0L + a2 * b;
			long double a = 1.0L - a1 * b;
			long double greb = d / c;

#if 0
			 printf("a1=%3d a2=%3d  a=%4.0Lf b=%4.0Lf c=%4.0Lf d=%4.0Lf"
				"  xlo=%7.5Lg  xhi=%7.5Lg  d/c=%Lg\n",
				na1, na2, a,b,c,d, xlo, xhi, greb);
#endif

			if ((1 == na1) && (2 == na2))
			{
				// Use alternative for the integral, since the other one
				// converges very poorly in this case. And v-v the other
				// one converges badly in all the other cases.
				// printf("special case \n");
				long double rat = a/c;
				long double complex thi, tlo;
				thi = special_eff(s, greb, xhi);
				tlo = special_eff(s, greb, xlo);
				sum += rat*(thi-tlo);

				rat = b/c;
				thi = special_eff(s-1.0L, greb, xhi);
				tlo = special_eff(s-1.0L, greb, xlo);
				sum += rat*(thi-tlo);
			}

			/* special case, where c=d=0, so general formula does not apply */
			else if (na1 == na2)
			{
				// printf("duude c=0 for a1=%d\n", na1);

				term = effo(s, xhi);
				sum += term;

				term = effo(s, xlo);
				sum -= term;
			}

			// General case
			else
			{
				long double rat = a/c;
				long double complex thi, tlo;
				thi = eff(s, greb, xhi);
				tlo = eff(s, greb, xlo);
				term = rat*(thi-tlo);
				// printf("thi=%g +I%g        tlo=%g+I%g\n\n", creal(thi), cimag(thi), creal(tlo), cimag(tlo));

				rat = b/c;
				thi = eff(s-1.0L, greb, xhi);
				tlo = eff(s-1.0L, greb, xlo);
				term += rat*(thi-tlo);
				// printf("thi=%g +I%g        tlo=%g+I%g\n\n", creal(thi), cimag(thi), creal(tlo), cimag(tlo));
				sum += term;

				if (cabs(term) < epsi)
				{
					cnt += cabs(term);
					break;
				}
				// printf("a1=%d a2=%d term=%g\n", na1, na2, cabs(term));
			}
		}
		// printf("a1=%d break at a2=%d\n", na1, na2);
	}
	*termerr = cnt;

	return sum;
}

long double complex zeta_12(long double complex s, double epsi, unsigned int a1max, double *errest)
{
	long double complex gral;

	gral = s/(s-1.0L);
	gral -= s* gral_s12(s, epsi, a1max, errest);

	return gral;
}

int main (int argc, char * argv[])
{
	int i;
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s <amax> <epsi>\n", argv[0]);
		exit(1);
	}
	int amax = atoi(argv[1]);
	double epsi = atof(argv[2]);

	long double complex ess = 0.5L;

	printf("#\n# max terms in summation=%d epsi=%g\n#\n", amax, epsi);
	for (i=0; i<500; i++)
	{
		double error_estimate;
		long double complex ans = zeta_12(ess, epsi, amax, &error_estimate);
		printf("%g	%12.10g	%12.10g	%g\n", cimag(ess), creal(ans), cimag(ans), error_estimate);
		fflush (stdout);
		ess += I*0.1L;
	}

	return 0;
}
