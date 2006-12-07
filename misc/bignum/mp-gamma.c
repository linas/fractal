
/*
 * mp-gamma.c
 *
 * Compute gamma function for various complex arguments
 *
 * Linas Vepstas December 2006
 */

#include <math.h>
#include <gmp.h>
#include "mp-binomial.h"
#include "mp-consts.h"
#include "mp-gamma.h"
#include "mp-misc.h"
#include "mp-trig.h"
#include "mp-zeta.h"

/* ================================================= */
/*
 * fp_lngamma -- compute log of gamma for real argument
 *
 * Uses simple, quickly converging algo-- A&S 6.1.33
 * Valid input must have -1 < x < 3
 */
static void reduced_lngamma (mpf_t gam, mpf_t ex, int prec)
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

/* ================================================= */

void fp_lngamma (mpf_t gam, const mpf_t z, int prec)
{
}

/* ================================================= */
/* 
 * gamma function, bug valid only for -1 < x < 3 
 */ 
static void reduced_gamma (mpf_t gam, mpf_t ex, int prec)
{
	reduced_lngamma (gam, ex, prec);
	fp_exp (gam, gam, prec);
}

/* ================================================= */
/* 
 * fp_gamma
 * Use pochhammer to get into range of 1 < z < 2,
 * Then use the reduced summation formula.
 */
void fp_gamma (mpf_t gam, const mpf_t z, int prec)
{
	mpf_t zee;
	mpf_init (zee);

	/* make a copy of the input arg NOW! */
	mpf_set (zee, z);
	
	/* double-presision used, this code doesn't need to 
	 * be all that accurate. */
	double flo = mpf_get_d (zee);
	if (flo > 2.0)
	{
		unsigned int intpart = (unsigned int) floor (flo-1.0);
		mpf_sub_ui (zee, zee, intpart);
		fp_poch_rising (gam, zee, intpart);
	}
	else if (flo < 1.0)
	{
		unsigned int intpart = (unsigned int) floor (2.0-flo);
		fp_poch_rising (gam, zee, intpart);
		mpf_ui_div (gam, 1, gam);

		mpf_add_ui (zee, zee, intpart);
	}
	else
	{
		mpf_set_ui (gam, 1);
	}

	mpf_t rgamma;
	mpf_init (rgamma);
	reduced_gamma (rgamma, zee, prec);
	
	mpf_mul (gam, gam, rgamma);

	mpf_clear (zee);
	mpf_clear (rgamma);
}

/* ==================  END OF FILE ===================== */
