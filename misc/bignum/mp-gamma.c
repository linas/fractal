
/*
 * mp-gamma.c
 *
 * Compute gamma function for various complex arguments
 *
 * Linas Vepstas December 2006
 */

#include <gmp.h>
#include "mp-consts.h"
#include "mp-gamma.h"
#include "mp-misc.h"
#include "mp-trig.h"
#include "mp-zeta.h"

/**
 * fp_lngamma -- compute log of gamma for real argument
 *
 * Uses simple, quickly converging algo-- A&S 6.1.33
 */
void fp_lngamma (mpf_t gam, mpf_t ex, int prec)
{
	int n;
	mpf_t z, zn, term;

	mpf_init (z);
	mpf_init (zn);
	mpf_init (term);

	/* make copy of input argument now! */
	mpf_set (z, ex);

	fp_log (gam, z, prec);
	mpf_neg (gam, gam);
	
	mpf_sub_ui (z, z, 1);
	mpf_mul (zn, z,z);

	/* Use 10^{-prec} for smallest term in sum */
	mpf_t maxterm;
	mpf_init (maxterm);
	fp_epsilon (maxterm, prec);
	
	n=2;
	while (1)
	{
		fp_zeta (term, n, prec);
		mpf_sub_ui (term, term, 1);
		mpf_mul (term, term, zn);
		mpf_div_ui (term, term, n);
		if (n%2)
		{
			mpf_sub (gam, gam, term);
		}
		else
		{
			mpf_add (gam, gam, term);
		}

		/* don't go no farther than this */
		mpf_abs(term, term);
		if (mpf_cmp (term, maxterm) < 0) break;

		mpf_mul (zn,zn,z);
		n++;
	}

	fp_euler_mascheroni (term, prec);
	mpf_sub_ui (term, term, 1);
	mpf_mul (term, term, z);
	mpf_sub (gam, gam, term);

	mpf_clear (z);
	mpf_clear (zn);
	mpf_clear (term);
}
