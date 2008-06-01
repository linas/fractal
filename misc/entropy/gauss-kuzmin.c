
/*
 * Gauss-Kuzmin entropy
 *
 * Linas Vepstas June 2008
 */

#include <math.h>
#include <stdio.h>

void entropy(void)
{
	long double p_tot, h;
	p_tot = 0.0;
	h = 0;
	int k;
	for (k=1; k<10000; k++)
	{
		long double p_k = k+1;
		p_k = 1.0 - 1.0/(p_k*p_k);
		p_k = -logl(p_k) / M_LN2;
		p_tot += p_k;
		h += p_k * logl(p_k);
		printf ("duude its %d %lg %lg %lg\n", k, p_k, p_tot, h);
	}
}

main ()
{
	entropy();
}
