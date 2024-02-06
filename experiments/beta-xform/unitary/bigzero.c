/*
 * bigzero.c
 *
 * Compute roots of polynomials.
 * Feb 2024
 */

#define _GNU_SOURCE
#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <anant/mp-zerofind.h>
#include <anant/mp-zeroiso.h>

#include "bigseries.c"

// -------------------------------------------------------

typedef struct {
	char* digs;
	int degree;
} bitsy;

static void isowrap(cpx_t f, int deriv, cpx_t z, void* args)
{
	bitsy* b = (bitsy*) args;
	ebz_deriv(f, deriv, z, b->digs, b->degree);
}

static void wrapper(cpx_t f, cpx_t z, int nprec, void* args)
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
 * Obtain roots of polynomial.
 * Roots are solutions of E = 0 = -1 + b_0 x + b_1 x^2 + b_k x^{k+1}
 * where degree == k+1
 * The bits are to be provided in the array "digs".
 * There are always "degree" roots. They are returned in the arrays
 * modulus and phase, which must have size greater than degree.
 * Returns number of roots actually found.
 */
int get_roots(char* digs, int degree, double* modulus, double* phase)
{
	bitsy b;
	b.digs = digs;
	b.degree = degree-1;  // off-by one because of definitions

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
	if (ndisk != degree)
		fprintf(stderr, "Error: expecting %d roots, got %d\n", degree, ndisk);

	mpf_t mod;
	mpf_init(mod);

	cpx_t zero;
	cpx_init(zero);
	cpx_t guess;
	cpx_init(guess);

#ifdef DEBUG
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
#endif

	for (int i=0; i<ndisk; i++)
	{
		cpx_set(guess, centers[i]);
		cpx_set_ui(e1, 1, 0);
		cpx_set_ui(e2, 0, 1);
		cpx_times_mpf(e1, e1, radii[i]);
		cpx_times_mpf(e2, e2, radii[i]);
		int rc = cpx_find_zero_r(zero, wrapper, guess, e1, e2, 25, 70, &b);
		if (0 != rc) printf("# Aieeeeeeeeeeeeeeeeeeee!\n");

		double re = cpx_get_re(zero);
		double im = cpx_get_im(zero);
		double r = sqrt(re*re + im*im);
		double phi = atan2(im, re) / M_PI;
		modulus[i] = r;
		phase[i] = phi;
	}

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

	return ndisk;
}

static int cmp(const void* vida, const void* vidb, void* args)
{
	double* angs = (double*) args;
	int ida = *((int*) vida);
	int idb = *((int*) vidb);
	double ma = angs[ida];
	double mb = angs[idb];
	return (ma < mb) ? -1: 1;
}

void sort_roots(int nroots, double* modulus, double* phase)
{
	// Sort zeros in angular order.
	int mpts[1000];
	for (int i=0; i< nroots; i++) mpts[i] = i;
	qsort_r(mpts, nroots, sizeof(int), cmp, phase);

	double tmp[1000];
	for (int i=0; i< nroots; i++) tmp[i] = modulus[mpts[i]];
	for (int i=0; i< nroots; i++) modulus[i] = tmp[i];

	for (int i=0; i< nroots; i++) tmp[i] = phase[mpts[i]];
	for (int i=0; i< nroots; i++) phase[i] = tmp[i];
}

// --------------------------- END OF FILE -----------------------
