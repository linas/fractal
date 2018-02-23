/*
 * fpcoef.c
 *
 * Numerical integration to obtain measure coefficients.
 * Plot dependence of coefficients w.r.t K
 *
 * Dec 2017
 * Feb 2018
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NOEMAIN
// #include "psi.c"
// #include "psibig.c"
#include "expan.c"

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s npts niter\n", argv[0]);
		exit(1);
	}

	int npts = atoi(argv[1]);
	int niter = atoi(argv[2]);

	printf("#\n# npts = %d depth=%d\n#\n", npts, niter);

	int nkv=800;
	for (int i=0; i< nkv; i++)
	{
		double un = ((double) i + 0.5) / ((double) nkv);
		// double K = 0.5 + 0.5*un;
		// double K = 0.62 - 0.12*un;
		double K = 1.0 - 0.441*un;

		find_midpoints(K, MAXN);
		// big_midpoints(K, 400, midpoints, MAXN);
		sequence_midpoints(K, MAXN);

#ifdef OLD_SCHOOL
		int siter = niter;
		if (K < 0.8) siter = 1.2*niter;
		if (K < 0.74) siter = 1.5*niter;
		if (K < 0.65) siter = 2*niter;
		if (K < 0.60) siter = 2.5*niter;
		if (K < 0.57) siter = 3*niter;
		if (K < 0.55) siter = 5*niter;
		if (K < 0.53) siter = 8*niter;
		if (K < 0.51) siter = 12*niter;
		double cf0 = alpha_n(K, 0, npts, siter);
		double cf1 = alpha_n(K, 1, npts, siter);
		double cf2 = alpha_n(K, 2, npts, siter);
		double cf3 = alpha_n(K, 3, npts, siter);
		printf ("%d %d	%g	%g	%g	%g	%g\n", i, siter, 2.0*K, cf0, cf1, cf2, cf3);
#endif
		int siter = niter;
		double cf = alpha_n(K, 1, npts, siter);
		siter++;
		double cfn = alpha_n(K, 1, npts, siter);
		while (1.0e-4 < fabs(cf-cfn))
		{
			cf = cfn;
			siter++;
			cfn = alpha_n(K, 1, npts, siter);
		}

		// Save - it will never get easier; only harder.
		niter = siter-1;

		double cf0 = alpha_n(K, 0, npts, siter);
		// double cf1 = alpha_n(K, 1, npts, siter);
		double cf1 = cfn;
		printf ("%d %d	%g	%g	%g", i, siter, 2.0*K, cf0, cf1);
		fflush(stdout);
		for (int ord=2; ord<11; ord++)
		{
			double cf = alpha_n(K, ord, npts, siter);
			printf("	%g", cf);
		}
		printf("\n");
		fflush(stdout);
	}
}
