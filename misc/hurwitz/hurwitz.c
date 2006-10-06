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

void hurwitz_zeta (double *hre, double *him, double sre, double sim, double q)
{
	int norder = 80;

	/* arrays storing values to be forward-differenced */
#define NARR 1000
	double refd[NARR];
	double imfd[NARR];
	double rens[NARR];
	double imns[NARR];

	int k;
	for (k=0; k<norder; k++)
	{
		double logkq = log(k+q);
		double mag = exp((1.0-sre) * logkq);
		refd[k] = mag * cos (sim*logkq);
		imfd[k] = mag * sin (sim*logkq);

		// printf ("its %d \t%g \t%g\n", k, refd[k], imfd[k]);
	}

	int n;
	for (n=0; n<norder; n++)
	{
		double rs=0.0;
		double is=0.0;
		double cyn = 1.0;
		for (k=0; k<=n; k++)
		{
			double bin = cyn*binomial (n,k);
			rs += bin * refd[k];
			is += bin * imfd[k];
			cyn = -cyn;
		}
		rens[n] = rs;
		imns[n] = is;
		printf ("its %d \t%g \t%g\n", n, rens[n], imns[n]);
	}
} 

main ()
{
	double q;
	double sre, sim;
	double hre, him;

	sre = 0.5;
	sim = 9.0;
	q = 0.3;
	hurwitz_zeta (&hre, &him, sre,sim, q);
}
