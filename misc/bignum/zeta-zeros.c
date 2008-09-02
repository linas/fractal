
/*
 * zeta-zeros.c
 *
 * Find and print zeta-function zeros by brute force.
 * Not very elegent, usful for low-end stuff.
 *
 * Linas Vepstas September 2008
 */

#include <gmp.h>

#include "mp-complex.h"
#include "mp-zeta.h"

void get_zeros()
{
	int prec = 30;

	mpf_t tee;
	mpf_init (tee);
	mpf_set_ui(tee, 14);

	cpx_t ess;
	cpx_init (ess);
	cpx_set (ess, 

	while(1)
	{
	}
	
	cpx_borwein_zeta(cpx_t zeta, const cpx_t ess, int prec);
}

main()
{
	get_zeros();
}
