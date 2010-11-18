/*
 * Eigenfunctions of the dyadic sawtooth, based on the 
 * Walsh functions. And integrals thereof.
 *
 * Linas Vepstas May 2010
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "mp-quest.h"

#define BITLEN 32

/* Stucture to hold bit-shift states, needed to compute eigenfunctions */
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
 * Count the total number of bits (zero or one) in the  integer n
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
	mpf_t ex;
	static int is_initialized=0;
	static mpf_t half;

	if (0 == is_initialized)
	{
		is_initialized = 1;
		mpf_init(half);
		mpf_set_ui(half, 1);
		mpf_div_ui(half, half, 2);
	}

	/* Compute the floor of x, and subtract it */
	mpf_init(ex);
	mpf_floor(ex, x);
	mpf_sub(result, x, ex);

	if (0 < mpf_cmp(half, result))
	{
		mpf_set_ui(result, 1);
		mpf_clear(ex);
		return;
	}
	
	mpf_set_si(result, -1);
	mpf_clear(ex);
}

/**
 * Triangle wave function function, returns integral of the 
 * square wave function (i.e. tent rising to 1/2 for x< 1/2 
 * and dropping to zero for x> 1/2
 * Repeats with period 1.
 */
void tent_1(mpf_t result, mpf_t x)
{
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
	mpf_floor(ex, x);
	mpf_sub(result, x, ex);

	if (0 < mpf_cmp(half, result))
	{
		mpf_clear(ex);
		return;
	}
	
	mpf_ui_sub(result, 1, result);
	mpf_clear(ex);
}

/**
 * Returns square wave of frequency 2^(n-1)
 */
void step_n(mpf_t result, mpf_t x, int n)
{
	mpf_t ex;

	if (1 == n)
	{
		step_1(result, x);
		return;
	}

	// Make copy of x, in case 'result' is same storage location as x.
	mpf_init(ex);
	mpf_mul_2exp(ex, x, n-1);
	step_1(result, ex);
	mpf_clear(ex);
}

/**
 * Returns integral of the square wave of frequency 2^(n-1)
 * The slope of the returned function is either +1 or -1.
 */
void igral_step_n(mpf_t result, mpf_t x, int n)
{
	if (1 == n)
	{
		tent_1(result, x);
		return;
	}

	mpf_mul_2exp(result, x, n-1); 
	tent_1(result, result);
	mpf_div_2exp(result, result, n-1);
}

/**
 * Implement the n'th walsh function.
 * n must be less than 2^32 or 2^64, depending on the
 * BITLEN setting.
 */
void walsh(mpf_t result, mpf_t x, unsigned long n)
{
	int i;
	mpf_t term, ex;
	mpf_init(term);
	mpf_init(ex);

	// Make copy of x, in case 'result' shares same storage location.
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

	mpf_clear(term);
	mpf_clear(ex);
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

	mpf_clear(term);
	mpf_clear(ex);
}

/**
 * Implement the n'th blancmange based on the n'th walsh function
 * w must be greater than zero.
 * n must be less that 2^32 or 2^64
 * prec is the desired number of decimal palces of precision.
 */
void blanc(mpf_t result, mpf_t w, mpf_t x, unsigned long n, int prec)
{
	double lw;
	int i, nterms;
	mpf_t ex, term, tn, wn;
	mpf_init(ex);
	mpf_init(term);
	mpf_init(wn);
	mpf_init(tn);

	mpf_set(ex, x);
	mpf_set_ui(wn, 1);
	mpf_set_ui(tn, 1);
	mpf_set_ui(result, 0);

	// Compute the required number of terms
	lw = - log(mpf_get_d(w));
	nterms = ((double) prec) * 2.3026 / lw;
	nterms += 4; // for good luck.

	for (i=1; i<=nterms; i++)
	{
		mpf_mul(term, tn, ex);
		walsh(term, term, n);
		mpf_mul(term, term, wn);
		mpf_add(result, result, term);

		mpf_mul(wn, wn, w);
		mpf_mul_ui(tn, tn, 2);
	}

	mpf_clear(ex);
	mpf_clear(term);
	mpf_clear(wn);
	mpf_clear(tn);

}

/**
 * Integral of the n'th blancmange based on the n'th walsh function.
 * By "blancmange" we mean the sum over powers of w i.e.
 * sum_k w^k f_n(2^k x)  where f_n is the n'th walsh func.
 * as usual, w should be less than 1 to strictly converge.
 * (converges conditinally for w less than 2).
 *
 * n must be less that 2^32 or 2^64
 * prec is the number of decimal digits of desired precision.
 */
void igral_blanc(mpf_t result, mpf_t w, mpf_t x, unsigned long n, int prec)
{
	double lw;
	int i, nterms;
	mpf_t term, ex, tn, wn;
	mpf_init(term);
	mpf_init(ex);
	mpf_init(wn);
	mpf_init(tn);

	mpf_set(ex, x);
	mpf_set_ui(wn, 1);
	mpf_set_ui(tn, 1);
	mpf_set_ui(result, 0);

	// Compute the required number of terms
	lw = - log(0.5 * mpf_get_d(w));
	nterms = ((double) prec) * 2.3026 / lw;
	nterms += 4; // for good luck.

	for (i=1; i<=nterms; i++)
	{
		/* Compute the integral of the n'th walsh function, for
		 * argument 2^(i-1)*x.  Remember to divide by 2^(i-1)*x
		 * afterwards, as otherwise integral is off by this factor.
		 * Altnerately, we could have multiplied n by 2, but this
		 * would have overflowed more easily, since n is limited.
		 */
		mpf_mul(term, tn, ex);
		igral_walsh(term, term, n);
		mpf_div(term, term, tn);    // Remember to normalize! Dohh!

		mpf_mul(term, term, wn);
		mpf_add(result, result, term);

		mpf_mul(wn, wn, w);
		mpf_mul_ui(tn, tn, 2);
	}

	mpf_clear(term);
	mpf_clear(ex);
	mpf_clear(wn);
	mpf_clear(tn);
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
#ifdef PRINT_WALSH_COEFFS_A1
// The below was used to graph the coefficient A1 for the gkw paper.
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

/*
 * Return a shift strcutre appropriate for this n, w. 
 *
 * The shifts are precomputed and cached, and can be re-used; 
 * WARNING, however, w MUST NOT CHAGE, else old, incorrect cached
 * values will be returned.
 */
Shifts * shift_cache(unsigned long n, mpf_t w)
{
	static unsigned int arraylen = 0;
	static Shifts **shifts = NULL;

	if ((n < arraylen) && shifts[n]) return shifts[n];

	/* Expand the array size if needed */
	if (arraylen <= n)
	{
		int more = n + n/8 + 1;
		shifts = realloc(shifts, more*sizeof (Shifts *));
		int i;
		for (i = arraylen; i < more; i++) shifts[i] = NULL;
		arraylen = more;
	}

	/* Alloc a new shift stuct and cache it */
	shifts[n] = (Shifts *) malloc(sizeof(Shifts));

	/* Initialize */
	get_shifts(shifts[n], n);

	double w_f = mpf_get_d(w);
	get_coeffs(shifts[n], w_f);

	return shifts[n];
}

/**
 * Compute the eigenfunction of the dyadic sawtooth, associated with w
 * prec is the desired decimal precision.
 */
void eigenfunc(mpf_t result, mpf_t w, mpf_t x, unsigned long n, int prec)
{
	int i;
	mpf_t ex, term, ak;
	mpf_init(ex);
	mpf_init(term);
	mpf_init(ak);

	mpf_set(ex, x);
	blanc(result, w, ex, n, prec);

	Shifts *sh = shift_cache(n, w);

	for (i=1; i<sh->bitlen; i++)
	{
		walsh(term, ex, sh->m_k[i]);
		mpf_set_d(ak, sh->a_k[i]);
		mpf_mul(term, term, ak);
		mpf_add(result, result, term);
	}

	mpf_clear(ex);
	mpf_clear(term);
	mpf_clear(ak);
}

/**
 * Compute the integral of the n'th eigenfunction of the dyadic sawtooth,
 * associated with w.  This is as given in the paper gkw.pdf with
 * the primary part given by the blancmange, plus a finite series
 * of "corrections" to turn this into an eigenfunc.
 *
 * Integral is from zero to x.
 *
 * prec is the number of decimal places of desired precision.
 */
void
igral_eigenfunc(mpf_t result, mpf_t w, mpf_t x,
                unsigned long n, int prec)
{
	int i;
	mpf_t ex, term, ak;
	mpf_init(ex);
	mpf_init(term);
	mpf_init(ak);

	mpf_set(ex, x);
	igral_blanc(result, w, ex, n, prec);

	Shifts *sh = shift_cache(n, w);

	for (i=1; i<sh->bitlen; i++)
	{
		igral_walsh(term, ex, sh->m_k[i]);
		mpf_set_d(ak, sh->a_k[i]);
		mpf_mul(term, term, ak);
		mpf_add(result, result, term);
	}

	mpf_clear(ex);
	mpf_clear(term);
	mpf_clear(ak);
}

/*
 * Compute an experimental linear combo of eigenfuncs, just to 
 * see what happens.  Under construction.
 */
void
sum_of_igral_eigenfunc(mpf_t result, mpf_t w, mpf_t x, int prec)
{
	unsigned int i, n;
	mpf_t ex, term, a_n;
	mpf_init(ex);
	mpf_init(term);
	mpf_init(a_n);

	mpf_set(ex, x);
	mpf_set_ui(result, 0);
	n = 1;

	for (i=1; i<4; i+=2)
	{
		n = i;
		igral_eigenfunc(term, w, ex, n, prec);

		int ln = bitlength(n);
		n = 1<<ln;

		mpf_set_ui(a_n, n);
		mpf_div_ui(a_n, a_n, ln);
		mpf_mul(term, term, a_n);

		mpf_add(result, result, term);

		// n = (1<<i) + 1;
	}

	mpf_clear(ex);
	mpf_clear(term);
	mpf_clear(a_n);
}

#ifdef PRINT_WALSH_COEFFS_A1
/**
 * This was used to generate the graph of a_1 for the the GKW paper.
 * the actual print statement is in the get_coeffs() function.
 * To get this to work, be sure to touch the #define in functioin above.
 */
int main (int argc, char * argv[])
{
	double  w_f;
	int n = 5;
	int prec, nbits;
	Shifts shifts;

	if (3 > argc)
	{
		fprintf(stderr, "Usage: %s <decimal-precision> <w>\n", argv[0]);
		exit(1);
	}

	/* prec is decimal-places of precision */
	prec = 50;
	prec = atoi(argv[1]);

	/* Set the precision (number of binary bits) */
	nbits = 3.3*prec;
	mpf_set_default_prec (nbits);

	/* Other misc args */
	w_f = atof(argv[2]);

	for (n=1; n<=200; n++)
	{
		get_shifts(&shifts, n);
		get_coeffs(&shifts, w_f);
	}
	return 0;
}
#endif

int main (int argc, char * argv[])
{
	mpf_t r, x, y, step, w;
	double r_f, x_f, y_f, f_f, w_f, yf_prev;
	int n = 5;
	int i, npts;
	int prec, nbits;

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

	mpf_init(x);
	mpf_init(y);
	mpf_init(r);
	mpf_init(step);
	mpf_init(w);
	mpf_set_d(w, w_f);

#if 1
	// Create graphs of integral for a single n at a time.
	npts = 1233;
	mpf_set_ui(step, 1);
	mpf_div_ui(step, step, npts);

	mpf_set_ui(r, 0);
	yf_prev = 0.0;
	for (i=0; i<npts; i++)
	{
		// step_1(y, x);
		// step_n(y, x, 3);

		r_f = mpf_get_d(r);  // dyadic place

		mpf_set(x, r);
		// Want integrals convoluted with question mark.
		question_mark(x, r, prec);  

		x_f = mpf_get_d(x);
		// walsh(y, x, n);
		// blanc(y, w, x, n, prec);
		// igral_walsh(y, x, n);
		// igral_blanc(y, w, x, n, prec);
		eigenfunc(y, w, r, n, prec);
		y_f = mpf_get_d(y);

		// igral_eigenfunc(y, w, x, n, prec);
		sum_of_igral_eigenfunc(y, w, x, prec);
		f_f = mpf_get_d(y);

		double delta = f_f - yf_prev;
		yf_prev = f_f;
		delta *= ((double) npts);

		// r == dyadic place
      // x == rational place
      // y == uhhh, crap, actually. Nonsense
      // f == inegral of eigenfunc in rational place
		// delta == hand-build derivative of integral == eigenfunc of gkw.
		printf("%d	%f	%f	%g	%g	%g\n", i, r_f, x_f, y_f, f_f, delta);

		mpf_add(r, r, step);
	}
#endif


#ifdef SHOW_DIVERGENCE_OF_DERIVATIVE
	/* The loop below explores a progressively narrower integration
	 * interval centered about x=2/3. i.e. an interaval from x to 
	 * x+ 2^(-n) which is exptect to .. whatever... the idea is graph and
	 * find out.
	 */
	mpf_set_ui(x, 2);
	mpf_div_ui(x, x, 3);
	npts=500;

	// graph of integral, 
	for (i=1; i<npts; i++)
	{
		int j;
		double prev = 1.0;
		double cur;
		printf("%d", i);
		for (j=0; j<6; j++)
		{
			mpf_set(y, x);
			igral_eigenfunc(y, w, y, i);
			f_f = mpf_get_d(y);

			mpf_set_ui(y, 1);
			mpf_div_2exp(y, y, n+j);
			mpf_add(y,y,x);
			igral_eigenfunc(y, w, y, i);
			f_f -= mpf_get_d(y);

#if 0
		cur = f_f;
		f_f /= prev;
		prev = cur;
#endif
			printf("	%g", f_f);
		}
		printf("\n");
	}
#endif

	mpf_clear(x);
	mpf_clear(y);
	mpf_clear(r);
	mpf_clear(step);
	mpf_clear(w);
	return 0;
}
