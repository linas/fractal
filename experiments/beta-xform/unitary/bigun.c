/*
 * bigun.c
 *
 * Compute "unitary" zeros using bignum.
 * Feb 2024
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <anant/mp-complex.h>
#include <anant/mp-zerofind.h>

/**
 * Compute binary digit sequence for a given beta value.
 *
 * digs: pointer to byte array in which to store bit sequence.
 * maxn: how many digits to compute.
 *
 * Each digit records whether 0.5<=midp (true) or not (false).
 * This is the Renyi-Parry sequence. The zeroth bit is always one.
 * the first bit is 0.5 <= (m_1 = T(beta/2))
 */
void
gen_bitseq(mpf_t beta, char* digs, int maxn)
{
	mpf_t half;
	mpf_init(half);
	mpf_set_ui(half, 1);
	mpf_div_ui(half, half, 2);

	mpf_t midp;
	mpf_init(midp);
	mpf_set(midp, beta);
	mpf_div_ui(midp, midp, 2);

	/* OK, now start iterating the map */
	for (int i=0; i < maxn; i++)
	{
		if (0 <= mpf_cmp(midp, half))
		{
			mpf_sub(midp, midp, half);
			digs[i] = 1;
		}
		else digs[i] = 0;
		mpf_mul(midp, midp, beta);
	}

	mpf_clear(half);
	mpf_clear(midp);
}

/**
 * Compute zeta polynomial. Roots of this poly will be the desired
 * zeros that we hope are accumulating.  This is the E(beta;z) thing.
 *
 * sum: returned value
 * zeta: location at which to evaluate.
 * digs: arry of digits to use
 * order: max order of the poly. There must be at least this many
 *     binary digits available.
 */
void ebz(cpx_t sum, cpx_t zeta, char* digs, int order)
{
	cpx_t zetan;
	cpx_init(zetan);
	cpx_set_ui(zetan, 1, 0);

	cpx_set_ui(sum, 0, 0);

	// Do the first k-1 of them.  The last one is always one.
	for (int i=0; i<order; i++)
	{
		if (digs[i]) cpx_add(sum, sum, zetan);
		cpx_mul(zetan, zetan, zeta);
	}
	// The final bit is always one.
	cpx_add(sum, sum, zetan);

	// times zeta
	cpx_mul(sum, sum, zeta);

	// subtract one.
	cpx_sub_ui(sum, sum, 1, 0);

	cpx_clear(zetan);
}

typedef struct {
	char* digs;
	int order;
} bitsy;

void wrapper(cpx_t f, cpx_t z, int nprec, void* args)
{
	bitsy* b = (bitsy*) args;
	ebz(f, z, b->digs, b->order);

	double re = cpx_get_re(z);
	double im = cpx_get_im(z);
	double r = sqrt(re*re + im*im);
	double phi = atan2(im, re) / (2.0*M_PI);
	double fre = cpx_get_re(f);
	double fim = cpx_get_im(f);
	double fmo = sqrt(fre*fre + fim*fim);
	printf("wrapni r=%f phi=%f fmodulus=%g\n", r, phi, fmo);
}

/**
 * Return zero of polynomial
 */
void survey(char* digs, int order)
{
	bitsy b;
	b.digs = digs;
	b.order = order;

	cpx_t e1;
	cpx_init(e1);
	cpx_t e2;
	cpx_init(e2);

	cpx_t zero;
	cpx_init(zero);
	cpx_t guess;
	cpx_init(guess);

	double scale = 6.4 / order;

	for (int i=0; i<2*order; i++)
	{
		double thet = 2.0* M_PI *i / ((double) 2*order);
		double co = cos(thet);
		double si = sin(thet);
		cpx_set_d(guess, 0.9*co, 0.9*si);

		cpx_set_d(e1, -0.5*scale*co, -0.5*scale*si);
		cpx_set_d(e2, scale*si, scale*co);

		int rc = cpx_find_zero_r(zero, wrapper, guess, e1, e2, 40, 80, &b);
		double re = cpx_get_re(zero);
		double im = cpx_get_im(zero);
		double r = sqrt(re*re + im*im);
		double phi = atan2(im, re) / (2.0*M_PI);
		printf("found one rc=%d r=%f phi=%f  x=%f  y=%f\n", rc, r, phi, re, im);
	}

	cpx_clear(e1);
	cpx_clear(e2);
	cpx_clear(guess);
	cpx_clear(zero);
}

int main(int argc, char* argv[])
{
	int bprec = 500;
	mpf_set_default_prec(bprec);

	// Set beta to exactly 1.6
	mpf_t beta;
	mpf_init(beta);
	mpf_set_ui(beta, 16);
	mpf_div_ui(beta, beta, 10);

#define NBITS 400
	char bitseq[NBITS];
	gen_bitseq(beta, bitseq, NBITS);
	for (int i=0; i<70; i++)
		printf("%d", bitseq[i]);
	printf("\n");

	survey(bitseq, 25);
}
