/*
 * Eigenfunctions of the dyadic sawtooth, based on the 
 * Walsh functions.
 *
 * Linas Vepstas May 2010
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#define BITLEN 32

typedef struct
{
	unsigned long m;
	int bitlen;
	unsigned long m_k[BITLEN];
	int sigma_k[BITLEN]; 
	double w;
	double a_k[BITLEN];
} Shifts;


/**
 * Square wave function function, returns +1 if x< 1/2 and returns -1 if x> 1/2
 * Repeats with period 1.
 */
void step_1(mpf_t result, mpf_t x)
{
	mpz_t intpart;
	mpf_t ex;
	static int is_init=0;
	static mpf_t half;

	if (0 == is_init)
	{
		is_init = 1;
		mpf_init(half);
		mpf_set_ui(half, 1);
		mpf_div_ui(half, half, 2);
	}

	mpf_init(ex);
	mpf_set(ex, x);
	mpz_init(intpart);
	mpz_set_f(intpart, ex);
	mpf_set_z(result, intpart);
	mpf_sub(result, ex, result);

	if (0 < mpf_cmp(half, result))
	{
		mpf_set_ui(result, 1);
		return;
	}
	
	mpf_set_si(result, -1);
}

/**
 * returns square wave of frequency 2^(n-1)
 * n must be less than 32/64
 * x must be in range of zero, one
 */
void step_n(mpf_t result, mpf_t x, int n)
{
	mpf_t ex;
	unsigned long tp;

	if (1 == n)
	{
		step_1(result, x);
		return;
	}

	mpf_init(ex);
	mpf_set(ex, x);

	tp = 1<<(n-1);
	mpf_set_ui(result, tp);
	mpf_mul(result, result, ex); 
	step_1(result, result);
}

/**
 * Implement the n'th walsh function
 * n must be less that 2^32 or 2^64
 */
void walsh(mpf_t result, mpf_t x, unsigned long n)
{
	int i;
	mpf_t term, ex;
	mpf_init(term);
	mpf_init(ex);
	mpf_set(ex, x);

	mpf_set_ui(result, 1);
	for (i=1; i<=BITLEN; i++)
	{
		if (n & 0x1)
		{
			step_n(term, ex, i);
			mpf_mul(result, result, term);
		}

		n >>= 1;
	}
}

/**
 * Implement the n'th blancmange based on the n'th walsh function
 * n must be less that 2^32 or 2^64
 */
void blanc(mpf_t result, mpf_t w, mpf_t x, unsigned long n)
{
	int i;
	mpf_t term, tn, wn;
	mpf_init(term);
	mpf_init(wn);
	mpf_init(tn);

	mpf_set_ui(wn, 1);
	mpf_set_ui(tn, 1);
	mpf_set_ui(result, 0);

	// XXX should be some variiable number of terms
	for (i=1; i<=60; i++)
	{
		mpf_mul(term, tn, x);
		walsh(term, term, n);
		mpf_mul(term, term, wn);
		mpf_add(result, result, term);

		mpf_mul(wn, wn, w);
		mpf_mul_ui(tn, tn, 2);
	}
}

static int bitcount(unsigned long n)
{
	int cnt = 0;
	while (n)
	{
		if (n & 0x1) cnt ++;
		n >>= 1;
	}
	return cnt;
}

void get_shifts(Shifts *sh, unsigned long n)
{
	int i, m;
	bzero(sh, sizeof(Shifts));
	sh->m = n;

	/* Count total number of bits in n */
	i = 0;
	m = n;
	while (m)
	{
		i++;
		m >>=1;
	}
	sh->bitlen = i;

	/* Store shift states */
	i = sh->bitlen;
	m = n;
	while (m)
	{
		sh->m_k[i] = m;
		i--;
		m >>=1;
	}

	/* Compute the sign bits */
	i = sh->bitlen;
	m = n;
	while (m)
	{
		int cnt = bitcount(m);
		if (cnt%2 == 0) sh->sigma_k[i] = 1;
		else sh->sigma_k[i] = -1;
		i--;
		m >>=1;
	}

printf ("duuude n=%lu bitlen=%d shifts=%lu %lu %lu %lu sigmas=%d %d %d %d\n", 
	n, sh->bitlen,
	sh->m_k[1], sh->m_k[2], sh->m_k[3], sh->m_k[4],
	sh->sigma_k[1], sh->sigma_k[2], sh->sigma_k[3], sh->sigma_k[4]);
}

int main (int argc, char * argv[])
{
	mpf_t x, y, step, w;
	double x_f, y_f, f_f, w_f;
	int n = 5;
	int i, npts;
	int prec, nbits;
	Shifts shifts;

	if (4 > argc)
	{
		fprintf(stderr, "Usage: %s <decimal-precision> <n> <w>\n", argv[0]);
		exit(1);
	}

	/* prec is decimal-places of precision */
	prec = 50;
	prec = atoi(argv[1]);

	/* Set the precision (number of binary bits) */
	nbits = 3.3*prec;
	mpf_set_default_prec (nbits);

	/* Other misc args */
	n = atoi(argv[2]);
	w_f = atof(argv[3]);

	get_shifts(&shifts, n);

	mpf_init(x);
	mpf_init(y);
	mpf_init(step);
	mpf_init(w);
	mpf_set_d(w, w_f);

	npts = 1233;
	mpf_set_ui(step, 1);
	mpf_div_ui(step, step, npts);

	mpf_set_ui(x, 0);
	for (i=0; i<npts; i++)
	{
		// step_1(y, x);
		// step_n(y, x, 3);

		x_f = mpf_get_d(x);
		walsh(y, x, n);
		y_f = mpf_get_d(y);
		blanc(y, w, x, n);
		f_f = mpf_get_d(y);

		printf("%d	%f	%g	%g\n", i, x_f, y_f, f_f);

		mpf_add(x, x, step);
	}

	return 0;
}
