/* 
 * hurwitz.c
 *
 * Compute the Hurwitz zeta function for arbitrary complex argument
 *
 * Linas Vepstas October 2006
 */

#include <math.h>
#include <stdio.h>

#include <gmp.h>
#include "mp-binomial.h"
#include "mp-cache.h"
#include "mp-complex.h"
#include "mp-misc.h"
#include "mp-trig.h"

/* =========================================================== */

/* Compute n'th forward difference of (k+q)^{s} */
static void forward_diff_diri (cpx_t fin, int n, mpf_t q, cpx_t ess, int prec)
{
	mpf_set_ui (fin[0].re, 0);
	mpf_set_ui (fin[0].im, 0);

	mpz_t ibin;
	mpz_init (ibin);
	
	mpf_t bin;
	mpf_init (bin);

	cpx_t diri;
	cpx_init (diri);

	int k;
	for (k=0; k<=n; k++)
	{
		i_binomial (ibin, n,k);
		mpf_set_z (bin, ibin);

		fp_pow_rc (diri, k, q, ess, prec);

		cpx_mul_mpf (diri, diri, bin);

		if (0 == k%2)
		{
			cpx_add (fin, fin, diri);
		}
		else
		{
			cpx_sub (fin, fin, diri);
		}
	}

	mpz_clear (ibin);
	mpf_clear (bin);
	cpx_clear (diri);
}

/* =========================================================== */
/* A brute-force summation using Hasse formula, 
 * for complex s, real q.
 *
 * Unfortunately, the convergence is unbearably slow, seems to be logarithmic....
 */

void hurwitz_zeta(cpx_t zeta, cpx_t ess, mpf_t q, int prec)
{
	int norder = prec;

	mpf_set_ui (zeta[0].re, 0);
	mpf_set_ui (zeta[0].im, 0);

	/* smo contains value of (1-s) */
	cpx_t smo;
	cpx_init (smo);
	cpx_neg (smo, ess);
	cpx_add_ui (smo, smo, 1, 0);

	/* os contains value of 1/(s-1) */
	cpx_t os;
	cpx_init (os);
	mpf_sub_ui (os[0].re, ess[0].re, 1);
	mpf_set (os[0].im, ess[0].im);
	cpx_recip (os, os);

	/* containes finite difference */
	cpx_t fd;
	cpx_init (fd);

	mpf_t on;
	mpf_init (on);

	int n;
	for (n=0; n<norder; n++)
	{
		forward_diff_diri (fd, n, q, smo, prec);
		
		mpf_set_ui (on, 1);
		mpf_div_ui (on, on, n+1);
		
		cpx_mul_mpf (fd, fd, on);
printf ("duude %d ", n);
fp_prt (" ", fd[0].re);
		cpx_add (zeta, zeta, fd);
cpx_mul (fd, zeta, os);
fp_prt ("   ", fd[0].re);
printf ("\n");
	}

	cpx_mul (zeta, zeta, os);

	cpx_clear (os);
	
	mpf_clear (on);
	cpx_clear (fd);
} 

int main ()
{
	int prec = 180;

	/* Set the precision (number of binary bits) */
	mpf_set_default_prec (3.3*prec+600);

	cpx_t ess, zeta;
	cpx_init (ess);
	cpx_init (zeta);

	cpx_set_d (ess, 0.5, 4.0);
	cpx_set_d (ess, 2.0, 0.0);
	
	mpf_t que;
	mpf_init (que);
	mpf_set_d (que,1.0);
	
	hurwitz_zeta (zeta, ess, que, prec);

	fp_prt ("its ", zeta[0].re);
	printf ("\n");

	return 0;
}
