/**
 * yarh.c
 * Port of swap.C to bignum
 *
 * FUNCTION:
 * Integral the permuation group of continued fractions
 * Expect to get Riemann zeta in the Gauss map case
 * and that is what we seem to get ... need high integration 
 * order though to get anything on the r=1/2 axis ... 
 *
 * Results written up in yarh.lyx
 *
 * Linas Feb 2005
 * Linas Vepstas December 2010
 */

#include <stdio.h>
#include <stdlib.h>

#include "mp-complex.h"
#include "mp-trig.h"

/**
 * Swap the first and second digits of the continued fraction
 *
 * nprec == number of deciman places of precision.
 */
void swap_1_2 (mpf_t y, mpf_t x, int nprec)
{
	static int init = 0;
	static mpf_t zero;
	if (!init)
	{
		init = 1;
		mpf_init(zero);
		mpf_set_ui(zero, 0);
	}

	/* a1 and a2 are the first two digitis of the 
	 * continued fraction */
	mpf_t ox, a1, a2;
	mpf_init (ox);
	mpf_init (a1);
	mpf_init (a2);

	mpf_ui_div (ox, 1, x);
	mpf_floor (a1, ox);
	mpf_sub(y, ox, a1);
	if (mpf_eq (y, zero, 3.32*nprec))
		goto done;

	/* Now get the second digit */
	mpf_ui_div (ox, 1, y);
	mpf_floor (a2, ox);
	mpf_sub(y, ox, a2);
	
	/* re-assemble the continued fraction */
	mpf_add(ox, y, a1);
	mpf_ui_div(y, 1, ox);
	mpf_add (ox, y, a2);
	mpf_ui_div(y, 1, ox);

done:
	mpf_clear (ox);
	mpf_clear (a1);
	mpf_clear (a2);
}

/**
 * The integrand, which is swap(x) * x^s 
 */
void integrand(cpx_t y, mpf_t x, cpx_t s, int nprec)
{
	mpf_t perm;
	mpf_init (perm);

	swap_1_2 (perm, x, nprec);

	cpx_mpf_pow (y, x, s, nprec);
	cpx_times_mpf (y, y, perm);

	mpf_clear (perm);
}

/**
 * Compute single integral of the integrand.
 * actually compute 
 * zeta = s/(s-1) - s \int_0^1 swap(x) x^{s-1} dx 
 *
 * Done using the simplest Newton sum possible.
 */
void integral(cpx_t y, unsigned int nsteps, cpx_t s, int nprec)
{
	mpf_t step;
	mpf_init (step);

	mpf_clear (step);
}

int main (int argc, char * argv[])
{

	return 0;
}
