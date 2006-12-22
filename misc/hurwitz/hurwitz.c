/* 
 * hurwitz.c
 *
 * Compute the Hurwitz zeta function for arbitrary complex argument
 *
 * Linas Vepstas October 2006
 */

#include <math.h>
#include <stdio.h>

#include "binomial.h"

/* A brute-force summation using Hasse formula, 
 * for complex s, real q.
 *
 * Unfortunately, the convergence is slow on the critical strip,
 * and double precision is not enough to do anything useful here.
 */

void hurwitz_zeta (long double *phre, long double *phim, long double sre, long double sim, long double q)
{
	int norder = 80;

	/* arrays storing values to be forward-differenced */
#define NARR 1000
	long double refd[NARR];
	long double imfd[NARR];

	int k;
	for (k=0; k<norder; k++)
	{
		long double logkq = logl(k+q);
		long double mag = expl((1.0L-sre) * logkq);
		refd[k] = mag * cosl (sim*logkq);
		imfd[k] = mag * sinl (-sim*logkq);

		// printf ("its %d \t%Lg \t%Lg\n", k, refd[k], imfd[k]);
	}

	long double hre = 0.0;
	long double him = 0.0;
	int n;
	for (n=0; n<norder; n++)
	{
		long double rs=0.0L;
		long double is=0.0L;
		long double cyn = 1.0L;
		for (k=0; k<=n; k++)
		{
			long double bin = cyn*binomial (n,k);
			rs += bin * refd[k];
			is += bin * imfd[k];
			cyn = -cyn;
		}

		long double on = 1.0L/(n+1.0L);
		hre += on * rs;
		him += on * is;
		printf ("its %d \t%Lg \thre=%Lg \t%Lg \thim=%Lg\n", n, rs, hre,is, him);
	}
} 

main ()
{
	long double q;
	long double sre, sim;
	long double hre, him;

	sre = 0.5;
	sim = 2.0;
	q = 0.3;

	sre = 2.0;
	sim = 0;
	q = 1.0;

	hurwitz_zeta (&hre, &him, sre,sim, q);
}
