/* 
 * hurwitz.c
 *
 * Compute the Hurwitz zeta function for arbitrary complex argument
 *
 * Linas Vepstas October 2006
 */

#include <math.h>
#include <stdio.h>

/* A brute-force summation using Hasse formula, 
 * for complex s, real q. */

void hurwitz_zeta (double *hre, double *him, double sre, double sim, double q)
{
	int order = 40;

	/* arrays storing values to be forward-differenced */
#define NARR 100
	double refd[NARR];
	double imfd[NARR];

	int k;
	for (k=0; k<norder; k++)
	{
		double logkq = log(k+q);
		double mag = exp((1.0-sre) * logkq);
		refd[k] = mag * cos (sim*logkq);
		imfd[k] = mag * sin (sim*logkq);

		printf ("its %d %g %g\n", k, refd[k], imfd[k]);
	}
} 

main ()
{
	double q;
	double sre, sim;
	double hre, him;

	sre = 0.5;
	sim = 23.0;
	q = 0.3;
	hurwitz_zeta (&hre, &him, sre,sim, q);
}
