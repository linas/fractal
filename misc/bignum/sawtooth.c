
/*
 * sawtooth.c
 *
 * evaluate second sawtooth transfer operator coeffs
 * Linas Feb 2010
 */

#include <gmp.h>
#include <stdio.h>

void eig(mpq_t result, int k, int n)
{
	int m;
	int tk, tn, sk, sn;
	int num, deno;
	mpq_t term, tmp;
	mpz_t bin;

	if (k == n)
	{
		mpq_set_ui(result,1,1);
		return;
	}

	mpq_init(term);
	mpq_init(tmp);
	mpz_init(bin);

	mpq_set_ui(result,0,1);
	for (m=k+1; m<=n; m++)
	{
		int tm;
		eig(term, m,n);
		tm = 1<<m;  // 2^m
		mpq_set_ui(tmp, tm, 2*tm - 1);
		
		mpq_mul(term, term, tmp);

		mpz_bin_uiui (bin, m, k);
		mpq_set_z(tmp, bin);
		mpq_mul(term, term, tmp);
	
		mpq_add(result, term, term);
	}

	tn = 1<<(n+1);
	tk = 1<<(k+1);
	sn = 1;
	if (n%2 == 1) sn = -1;
	sk = 1;
	if (k%2 == 1) sk = -1;

	num = 2*sn*(tn-1)*(tk-1);
	deno = tk*(sk*(tk-1) - sn*(tn-1));
	// printf("duude k=%d num=%d den=%d\n", k, num, deno);
	mpq_set_si(tmp, num, deno);
	mpq_mul(result, result, tmp);

	mpq_canonicalize(result);
	
	mpq_clear(term);
	mpq_clear(tmp);
}

int
main(int argc, char * argv[])
{
	int k,n;
	mpq_t e;
	mpq_init(e);
	n = 3;

	for (k=0; k<=n; k++)
	{
		eig(e,k,n);
		printf("duude %d %d = ", k,n);
		mpq_out_str(stdout, 10, e);
		printf("\n");
	}

	return 0;
}
