/*
 * bracket-steps.c
 *
 * Look at fronts as a sum over square waves.
 * Sounded like a good idea, but isn't going to work cause of
 * unpredictable leaders.
 *
 * Similar to front-dist.c
 *
 * December 2023
 */

#include "selfie.c"
#include "selfie-util.c"
#include "selfie-tree.c"

int main(int argc, char* argv[])
{
	// Print front values as plateaus.
	if (2 != argc) {
		fprintf(stderr, "Usage: %s <maxord>\n", argv[0]);
		exit(1);
	}

	unsigned long prev = 0;
	double sklast = 0.0;
	double fact = 1.0;

	// int nu = k+2
	int maxord = atoi(argv[1]);
	for (int k=0; k<maxord; k++)
	{
		fact *= k+2;
		unsigned long mstart = 1UL << k;
		unsigned long mend = 1UL << (k+1);
		for (unsigned long m=mstart; m<mend; m++)
		{
			double x = (((double) m-mstart) + 0.5) / ((double) mstart);
			long idx = good_index_map(m);
			double nrm = ((double) mstart);
			// nrm = nrm*nrm;
			nrm  = fact;
			double sk = ((double) idx) / nrm;

			printf("%ld	%ld	%g	%ld	%g\n", m, m-mstart, x, prev, sklast);
			printf("%ld	%ld	%g	%ld	%g\n", m, m-mstart, x, idx, sk);
			prev = idx;
			sklast = sk;
		}
		printf("%ld	%ld	%g	%ld	%g\n", mend, mend-mstart, 1.0, prev, sklast);
		printf("\n");
		printf("%ld	%ld	%g	%ld	%g\n", mend, mend-mstart, 0.0, prev, sklast);
		fflush(stdout);
	}
}
