/*
 * bigun.c
 *
 * Compute "unitary" zeros using bignum.
 * Feb 2024
 */

#include "bigzero.c"

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// -------------------------------------------------------

int main(int argc, char* argv[])
{
	int bprec = 500;
	mpf_set_default_prec(bprec);
	printf("#\n# Default prec=%d bits\n#\n", bprec);

	// Set beta to exactly 1.6
	mpf_t beta;
	mpf_init(beta);
	mpf_set_ui(beta, 16);
	mpf_div_ui(beta, beta, 10);

#define NBITS 400
	char bitseq[NBITS];
	gen_bitseq(beta, bitseq, NBITS);
	printf("#\n# ");
	for (int i=0; i<70; i++)
		printf("%d", bitseq[i]);
	printf("\n#\n");

	int maxd = 101;
	double* modulus = (double*) malloc(maxd*sizeof(double));
	double* phase = (double*) malloc(maxd*sizeof(double));
	for (int degree=2; degree<maxd; degree ++)
	{
		get_roots(bitseq, degree, modulus, phase);
		sort_roots(degree, modulus, phase);

		for (int i=0; i<degree; i++)
		{
			double r = modulus[i];
			double phi = phase[i];
			printf("%d	%f	%f\n", degree, r, phi);
		}
		printf("\n");
		fflush(stdout);
	}
}
