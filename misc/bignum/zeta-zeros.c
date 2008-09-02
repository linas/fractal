
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

	mpf_t tee, half;
	mpf_init (tee);
	mpf_init (half);

	mpf_set_ui(half, 1);
	mpf_div_ui(half, half, 2);

	mpf_set_ui(tee, 14);

	cpx_t zeta, ess;
	cpx_init (ess);
	cpx_init (zeta);
	cpx_set_mpf (ess, half, tee);

	while(1)
	{
		cpx_borwein_zeta(zeta, ess, prec);

		mpf_add (ess[0].im, ess[0].im, half);

		fp_prt ("its ", zeta[0].re);
		printf ("\n");
	}
	
}

int main(int argc, char * argv[])
{
	get_zeros();
	return 0;
}
