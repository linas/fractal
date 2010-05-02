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
 * Count the total number of bits in the  integer n
 * (ignoring the infinite string of zeros padding to the left)
 */
static int bitlength(unsigned long n)
{
	int cnt = 0;
	while (n)
	{
		cnt ++;
		n >>= 1;
	}
	return cnt;
}

/**
 * Count the number of 1 bits in the  integer n
 */
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

	/* Compute the floor of x, and subtract it */
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
 * Triangle wave function function, returns integral of the 
 * square wave function (i.e. tent rising to 1/2 for x< 1/2 
 * and dropping to zero for x> 1/2
 * Repeats with period 1.
 */
void tent_1(mpf_t result, mpf_t x)
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

	/* Compute the floor of x, and subtract it */
	mpf_init(ex);
	mpf_set(ex, x);
	mpz_init(intpart);
	mpz_set_f(intpart, ex);
	mpf_set_z(result, intpart);
	mpf_sub(result, ex, result);

	if (0 < mpf_cmp(half, result))
	{
		return;
	}
	
	mpf_ui_sub(result, 1, result);
}

/**
 * returns square wave of frequency 2^(n-1)
 * n must be less than 32/64
 * x must be positive
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
 * Returns integral of the square wave of frequency 2^(n-1)
 * The slope of the returned function is either +1 or -1.
 * n must be less than 32/64
 * x must be positive
 */
void igral_step_n(mpf_t result, mpf_t x, int n)
{
	mpf_t ex;
	unsigned long tp;

	if (1 == n)
	{
		tent_1(result, x);
		return;
	}

	mpf_init(ex);
	mpf_set(ex, x);

	tp = 1<<(n-1);
	mpf_set_ui(result, tp);
	mpf_mul(result, result, ex); 
	tent_1(result, result);
	mpf_div_ui(result, result, tp);
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
 * Implement the integral of the n'th walsh function
 * n must be less that 2^32 or 2^64
 */
void igral_walsh(mpf_t result, mpf_t x, unsigned long n)
{
	int i, bitlen;
	mpf_t term, ex;
	mpf_init(term);
	mpf_init(ex);
	mpf_set(ex, x);

	bitlen = bitlength(n);

	mpf_set_ui(result, 1);
	for (i=1; i<bitlen; i++)
	{
		if (n & 0x1)
		{
			step_n(term, ex, i);
			mpf_mul(result, result, term);
		}
		n >>= 1;
	}

	igral_step_n(term, ex, bitlen);
	mpf_mul(result, result, term);
}

/**
 * Implement the n'th blancmange based on the n'th walsh function
 * n must be less that 2^32 or 2^64
 */
void blanc(mpf_t result, mpf_t w, mpf_t x, unsigned long n)
{
	int i;
	mpf_t ex, term, tn, wn;
	mpf_init(ex);
	mpf_init(term);
	mpf_init(wn);
	mpf_init(tn);

	mpf_set(ex, x);
	mpf_set_ui(wn, 1);
	mpf_set_ui(tn, 1);
	mpf_set_ui(result, 0);

	// XXX should be some variable number of terms
	for (i=1; i<=60; i++)
	{
		mpf_mul(term, tn, ex);
		walsh(term, term, n);
		mpf_mul(term, term, wn);
		mpf_add(result, result, term);

		mpf_mul(wn, wn, w);
		mpf_mul_ui(tn, tn, 2);
	}
}

/**
 * Implement the integral of the n'th blancmange based on the n'th walsh function
 * n must be less that 2^32 or 2^64
 */
void igral_blanc(mpf_t result, mpf_t w, mpf_t x, unsigned long n)
{
	int i;
	mpf_t term, ex, tn, wn;
	mpf_init(term);
	mpf_init(ex);
	mpf_init(wn);
	mpf_init(tn);

	mpf_set(ex, x);
	mpf_set_ui(wn, 1);
	mpf_set_ui(tn, 1);
	mpf_set_ui(result, 0);

	// XXX should be some variable number of terms
	for (i=1; i<=60; i++)
	{
		mpf_mul(term, tn, ex);
		igral_walsh(term, term, n);
		mpf_mul(term, term, wn);
		mpf_add(result, result, term);

		mpf_mul(wn, wn, w);
		mpf_mul_ui(tn, tn, 2);
	}
}

/**
 * Compute the assorted shift states needed for constructing
 * the dyadic sawtooth eigenfunctions.
 */
void get_shifts(Shifts *sh, unsigned long n)
{
	int i, m;
	bzero(sh, sizeof(Shifts));
	sh->m = n;

	/* Count total number of bits in n */
	sh->bitlen = bitlength(n);

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

#if 0
// A debugging printf
printf ("duuude n=%lu bitlen=%d shifts=%lu %lu %lu %lu sigmas=%d %d %d %d\n", 
	n, sh->bitlen,
	sh->m_k[1], sh->m_k[2], sh->m_k[3], sh->m_k[4],
	sh->sigma_k[1], sh->sigma_k[2], sh->sigma_k[3], sh->sigma_k[4]);
#endif
}

/**
 * Compute the explicit coeffs needed for the dyadic sawtooth 
 * eigenfuncs
 */
void get_coeffs(Shifts *sh, double w)
{
	int n, j, k;
	double tk;

	sh->w = w;
	n = sh->bitlen;
	k = 1;
	sh->a_k[n-k] = 1.0 / w;

	/* Fill in the rest of them */
	tk = 2.0;
	for (k=2; k<n; k++)
	{
		double acc = 0.0;
		double tn = 1.0;

		/* Make the j sum */
		for (j=1; j<k; j++)
		{
			acc += tn * sh->a_k[n-j];
			tn *= 2.0;
		}

		acc *= 2.0 - w;
		acc += 1.0;
		acc *= sh->sigma_k[n] * sh->sigma_k[n-k+1];
		acc /= tk * w;

		sh->a_k[n-k] = acc;
		tk *= 2.0;
	}
#if 0
// The below was used to graph the coefficients for the gkw paper.
// printf ("duuude m=%lu bitlen=%d aks=%g %g %g %g\n", 
printf ("%lu	%d	%g	%g	%g	%g\n", 
	sh->m-1,
	sh->bitlen,
	sh->a_k[1], sh->a_k[2], sh->a_k[3], sh->a_k[4]);
printf ("%lu	%d	%g	%g	%g	%g\n", 
	sh->m,
	sh->bitlen,
	sh->a_k[1], sh->a_k[2], sh->a_k[3], sh->a_k[4]);
#endif
}

/**
 * Compute the eigenfunction of the dyadic sawtooth, associated with w
 */
void eigenfunc(mpf_t result, mpf_t w, Shifts *sh, mpf_t x, unsigned long n)
{
	int i;
	mpf_t ex, term, ak;
	mpf_init(ex);
	mpf_init(term);
	mpf_init(ak);

	mpf_set(ex, x);
	blanc(result, w, x, n);

	for (i=1; i<sh->bitlen; i++)
	{
		walsh(term, x, sh->m_k[i]);
		mpf_set_d(ak, sh->a_k[i]);
		mpf_mul(term, term, ak);
		mpf_add(result, result, term);
	}
}

/**
 * Compute the integral of the eigenfunction of the dyadic sawtooth, associated with w
 */
void igral_eigenfunc(mpf_t result, mpf_t w, Shifts *sh, mpf_t x, unsigned long n)
{
	int i;
	mpf_t ex, term, ak;
	mpf_init(ex);
	mpf_init(term);
	mpf_init(ak);

	mpf_set(ex, x);
	igral_blanc(result, w, x, n);

	for (i=1; i<sh->bitlen; i++)
	{
		igral_walsh(term, x, sh->m_k[i]);
		mpf_set_d(ak, sh->a_k[i]);
		mpf_mul(term, term, ak);
		mpf_add(result, result, term);
	}
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
	get_coeffs(&shifts, w_f);

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
		// blanc(y, w, x, n);
		// igral_walsh(y, x, n);
		// igral_blanc(y, w, x, n);
		// eigenfunc(y, w, &shifts, x, n);
		igral_eigenfunc(y, w, &shifts, x, n);
		f_f = mpf_get_d(y);

		printf("%d	%f	%g	%g\n", i, x_f, y_f, f_f);

		mpf_add(x, x, step);
	}

	return 0;
}
