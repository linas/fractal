/*
 * complex-root.c
 *
 * Compute complex roots.
 * Feb 2024
 */

#include "bigzero.c"
#include "selfie.c"

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
		double r = modulus[i];
		double phi = phase[i];
		printf("r= %f	phi= %f\n", r, phi);
	}
}
