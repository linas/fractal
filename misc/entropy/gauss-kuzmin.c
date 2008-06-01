
/*
 * Gauss-Kuzmin entropy
 *
 * Linas Vepstas June 2008
 */

#include <math.h>
#include <stdio.h>

#ifndef M_LN2l
#define M_LN2l      0.6931471805599453094172321214581766L  /* log_e 2 */
#endif


void entropy(void)
{
	long double p_tot, h;

	p_tot = 0.0;
	h = 0;
	unsigned long k;
	for (k=1; k<2023012000; k++)
	{
		long double p_k = k+1;
		p_k = 1.0L - 1.0L/(p_k*p_k);
		p_k = -logl(p_k) / M_LN2l;
		p_tot += p_k;
		h -= p_k * logl(p_k);
		if (k%112001 == 0) printf ("duude its %d %Lg \t %Lg \t %18.15Lg\n", k, p_k, 1.0L - p_tot, h);
	}
}

main ()
{
	entropy();
}
