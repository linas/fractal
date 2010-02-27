
/*
 * sawtooth.c
 *
 * evaluate second sawtooth transfer operator coeffs
 * Linas Feb 2010
 */

#include <gmp.h>

void eig(mpq_t result, int k, int n)
{
	int m;
	mpq_t term, tmp;
	mpz_t bin, one;

	if (k == n)
	{
		mpq_set_ui(result,1,1);
		return;
	}

	mpq_init(term);
	mpq_init(tmp);
	mpz_init(bin);
	mpz_init(one);
	mpz_set_ui(one, 1);

	mpq_set_ui(result,0,1);
	for (m=k+1; m<=n; m++)
	{
		int tm;
		eig(term, m,n);
		tm = 1<<m;  // 2^m
		mpq_set_ui(tmp, tm, 2*tm - 1);
		
		mpq_mul(term, term, tmp);

		mpz_bin_uiui (bin, m, k);
		mpq_set_z(tmp, bin, one);
		mpq_mul(term, term, tmp);
	
		mpq_add(result, term, term);
	}

	mpq_canonicalize(result);
	
	mpq_clear(term);
	mpq_clear(tmp);
}
