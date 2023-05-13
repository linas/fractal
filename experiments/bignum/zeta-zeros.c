
/*
 * zeta-zeros.c
 *
 * Find and print zeta-function zeros by brute force.
 * Not very elegent, usful for low-end stuff.
 *
 * Linas Vepstas September 2008
 */

#include <gmp.h>
#include <stdio.h>

#include "mp-complex.h"
#include "mp-misc.h"
#include "mp-zeta.h"

void get_zeros()
{
	int prec = 30;

	mpf_t tee, half, step;
	mpf_init (tee);
	mpf_init (step);
	mpf_init (half);

	mpf_set_ui(half, 1);
	mpf_div_ui(half, half, 2);

	mpf_set_ui(step, 1);
	mpf_div_ui(step, step, 10);

	mpf_set_ui(tee, 13);

	cpx_t zeta, ess;
	cpx_init (ess);
	cpx_init (zeta);
	cpx_set_mpf (ess, half, tee);

	while(1)
	{
		cpx_borwein_zeta(zeta, ess, prec);

		mpf_add (ess[0].im, ess[0].im, step);

		fp_prt ("its ", ess[0].im);
		fp_prt ("\t", zeta[0].re);
		fp_prt ("\t", zeta[0].im);
		printf ("\n");
	}
	
}

int main(int argc, char * argv[])
{
	get_zeros();
	return 0;
}
