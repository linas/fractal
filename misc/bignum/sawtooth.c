
/*
 * sawtooth.c
 *
 * evaluate second sawtooth transfer operator coeffs
 * Linas Feb 2010
 */

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

void eig(mpq_t result, int k, int n)
{
	int m;
	int tk, tn, sk, sn;
	int num, deno;
	mpq_t acc, term, tmp;
	mpz_t bin;

	if (k == n)
	{
		mpq_set_ui(result,1,1);
		return;
	}

	mpq_init(term);
	mpq_init(tmp);
	mpq_init(acc);
	mpz_init(bin);

	mpq_set_ui(acc,0,1);
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
	
		mpq_add(acc, acc, term);
	}
#if 0
	mpq_canonicalize(acc);
	printf("accm(%d %d) = ", k,n);
	mpq_out_str(stdout, 10, acc);
	printf("\n");
#endif

	tn = 1<<(n+1);
	tk = 1<<(k+1);
	sn = 1;
	if (n%2 == 1) sn = -1;
	sk = 1;
	if (k%2 == 1) sk = -1;

	num = 2*sn*(tn-1)*(tk-1);
	deno = tk*(sk*(tk-1) - sn*(tn-1));
	if (deno < 0) { deno = -deno; num=-num; }
	// printf("duude k=%d num=%d den=%d\n", k, num, deno);
	mpq_set_si(tmp, num, deno);
	mpq_mul(acc, acc, tmp);

	mpq_canonicalize(acc);

#if 0
	printf("eig(%d %d) = ", k,n);
	mpq_out_str(stdout, 10, acc);
	printf("\n");
#endif

	mpq_set(result, acc);

	mpq_clear(term);
	mpq_clear(tmp);
	mpq_clear(acc);
}

int
main(int argc, char * argv[])
{
	int k,n;
	mpq_t e;
	mpq_init(e);
	n = 3;

	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <n>\n", argv[0]);
		exit(1);
	}

	n = atoi(argv[1]);

	for (k=0; k<=n; k++)
	{
		eig(e,k,n);
		printf("eigenvec(%d %d) = ", k,n);
		mpq_out_str(stdout, 10, e);
		printf("\n");
	}

	return 0;
}
