/*
 * complex-root.c
 *
 * Compute complex roots of golden polynomial of a given index.
 * Feb 2024
 */

#include "bigzero.c"
#include "selfie.c"
#include "metrics.c"

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// -------------------------------------------------------

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s <index>\n", argv[0]);
		exit(1);
	}
	long idx = atol(argv[1]);

	int bprec = 100;
	mpf_set_default_prec(bprec);
	printf("Default prec=%d bits\n", bprec);

	unsigned long tno = 2UL * idx + 1UL;
	int order = bitlen(tno);

#define NBITS 63
	char bitseq[NBITS];
	for (int i=0; i<order; i++)
	{
		bitseq[i] = (tno & 1UL);
		tno >>= 1;
	}

	printf("Order=%d bitseq= ", order);
	for (int i=0; i<order; i++)
		printf("%d", bitseq[i]);
	printf("\n");

	double modulus[NBITS];
	double phase[NBITS];
	get_roots(bitseq, order, modulus, phase);
	sort_roots(order, modulus, phase);

	for (int i=0; i<order; i++)
	{
		double zr = modulus[i];
		double zphi = phase[i];
		double r = 1.0/zr;
		printf("r=%f	1/r= %f	phi= %f\n", r, zr, zphi);
	}
	printf("\n");

	double beta = golden_beta(idx);
	printf("beta = %f\n", beta);
	printf("\n");

	for (int i=1; i<order; i++)
	{
		double ur = modulus[i-1];
		double uphi = phase[i-1];
		double vr = modulus[i];
		double vphi = phase[i];

		ur /= beta;
		vr /= beta;
		double delt = poincare_iso(ur, uphi, vr, vphi);
		double kdelt = klein_iso(ur, uphi, vr, vphi);
		printf("i=%d delt=%f klein=%f\n", i, delt, kdelt);
	}
}
