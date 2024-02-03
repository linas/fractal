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
#include <anant/mp-zeroiso.h>

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

// -------------------------------------------------------
/**
 * Compute zeta polynomial. Roots of this poly will be the desired
 * zeros that we hope are accumulating.  This is the E(beta;z) thing.
 *
 * sum: returned value
 * zeta: location at which to evaluate.
 * digs: arry of digits to use
 * order: max order of the poly. There must be at least this many
 *     binary digits available. Degree of polynomial is one more than
 *     this. Sorry for the off-by-one, but so it goes.
 */
void ebz(cpx_t sum, cpx_t zeta, char* digs, int order)
{
	cpx_t zetan;
	cpx_init(zetan);
	cpx_set(zetan, zeta);

	cpx_set_ui(sum, 0, 0);

	// Do the first k-1 bits.  The last bit is always one.
	for (int i=0; i<order; i++)
	{
		if (digs[i]) cpx_add(sum, sum, zetan);
		cpx_mul(zetan, zetan, zeta);
	}
	// The final bit is always one.
	cpx_add(sum, sum, zetan);

	// subtract one.
	cpx_sub_ui(sum, sum, 1, 0);

	cpx_clear(zetan);
}

// -------------------------------------------------------
/**
 * Same as above, but computes the m'th derivative.
 * Needed as a test function for zero finding.
 */
void ebz_deriv(cpx_t sum, int mderiv, cpx_t zeta, char* digs, int order)
{
	cpx_set_ui(sum, 0, 0);
	if (order+1 < mderiv) return;

	mpf_t fact;
	mpf_init(fact);
	mpf_set_ui(fact, 1);

	cpx_t zetan;
	cpx_init(zetan);

	if (0 == mderiv)
		cpx_set(zetan, zeta);
	else
		cpx_set_ui(zetan, 1, 0);

	cpx_t term;
	cpx_init(term);

	// Leading factorial
	for (int i=1; i<mderiv; i++)
		mpf_mul_ui(fact, fact, i+1);

	// Do the first k-1 bits.  The last bit is always one.
	int mstart = mderiv-1;
	if (mstart < 0) mstart = 0;
	for (int i=mstart; i<order; i++)
	{
		if (digs[i])
		{
			cpx_times_mpf(term, zetan, fact);
			cpx_add(sum, sum, term);
		}
		cpx_mul(zetan, zetan, zeta);
		if (0 < mderiv)
		{
			mpf_mul_ui(fact, fact, i+2);
			mpf_div_ui(fact, fact, i-mderiv+2);
		}
	}
	// The final bit is always one.
	cpx_times_mpf(term, zetan, fact);
	cpx_add(sum, sum, term);

	// Subtract one.
	if (0 == mderiv)
		cpx_sub_ui(sum, sum, 1, 0);

	cpx_clear(zetan);
	cpx_clear(term);
	mpf_clear(fact);
}

// -------------------------------------------------------

typedef struct {
	char* digs;
	int degree;
} bitsy;

void isowrap(cpx_t f, int deriv, cpx_t z, void* args)
{
	bitsy* b = (bitsy*) args;
	ebz_deriv(f, deriv, z, b->digs, b->degree);
}

void wrapper(cpx_t f, cpx_t z, int nprec, void* args)
{
	bitsy* b = (bitsy*) args;
	ebz(f, z, b->digs, b->degree);

#if 0
	double re = cpx_get_re(z);
	double im = cpx_get_im(z);
	double r = sqrt(re*re + im*im);
	double phi = atan2(im, re) / (2.0*M_PI);
	double fre = cpx_get_re(f);
	double fim = cpx_get_im(f);
	double fmo = sqrt(fre*fre + fim*fim);
	printf("wrapni r=%f phi=%f fmodulus=%g\n", r, phi, fmo);
#endif
}

/**
 * Return zero of polynomial
 */
void survey(char* digs, int degree)
{
	bitsy b;
	b.digs = digs;
	b.degree = degree-1;  // off-by one because of defintions

	// Allocate disks.
	int narr = degree+5;
	cpx_t* centers = (cpx_t*) malloc(narr*sizeof(cpx_t));
	mpf_t* radii = (mpf_t*) malloc(narr*sizeof(mpf_t));
	for (int i=0; i< narr; i++)
	{
		cpx_init(centers[i]);
		mpf_init(radii[i]);
	}

	// Bounding box
	cpx_t e1;
	cpx_init(e1);
	cpx_t e2;
	cpx_init(e2);
	cpx_set_ui(e1, 2, 2);
	cpx_neg(e2, e1);

	int ndisk = cpx_isolate_roots(isowrap, degree, e2, e1, centers, radii, &b);

	mpf_t mod;
	mpf_init(mod);

	cpx_t zero;
	cpx_init(zero);
	cpx_t guess;
	cpx_init(guess);

	printf("Degree %d found %d disks\n", degree, ndisk);
	for (int i=0; i<ndisk; i++)
	{
		cpx_abs(mod, centers[i]);
		printf("Found disk %2d center= %f %f radius= %f mod=%f\n", i,
			cpx_get_re(centers[i]),
			cpx_get_im(centers[i]),
			mpf_get_d(radii[i]),
			mpf_get_d(mod));

		cpx_set(guess, centers[i]);
		cpx_set_ui(e1, 1, 0);
		cpx_set_ui(e2, 0, 1);
		cpx_times_mpf(e1, e1, radii[i]);
		cpx_times_mpf(e2, e2, radii[i]);
		int rc = cpx_find_zero_r(zero, wrapper, guess, e1, e2, 25, 70, &b);
		cpx_abs(mod, zero);
		printf("            zero rc=%d %f %f                  mod=%f\n", rc,
			cpx_get_re(zero),
			cpx_get_im(zero),
			mpf_get_d(mod));
		if (0 != rc) printf("Aieeeeeeeeeeeeeeeeeeee!\n");
	}
	printf("----\n");
	fflush(stdout);

	cpx_clear(e1);
	cpx_clear(e2);
	mpf_clear(mod);
	cpx_clear(guess);
	cpx_clear(zero);

	for (int i=0; i< narr; i++)
	{
		cpx_clear(centers[i]);
		mpf_clear(radii[i]);
	}
	free(centers);
	free(radii);
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

	for (int degree=3; degree<25; degree++)
	{
		survey(bitseq, degree);
	}
}
