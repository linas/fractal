/*
 * psimid.c
 * Print values of midpoints.
 *
 * February 2018
 */
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define NOMAIN
#include "psi.c"
#include "psibig.c"
#undef NOMAIN

double bisect(double Klo, double Khi, int which)
{
	double Kmi = 0.5 * (Klo + Khi);
	if (Khi - Klo < 1.0e-16) return Kmi;

	big_midpoints(Klo, 400, midpoints, which+4);
	sequence_midpoints(Klo, which+4);
	double mlo = midpoints[which];

	big_midpoints(Khi, 400, midpoints, which+4);
	sequence_midpoints(Khi, which+4);
	double mhi = midpoints[which];

	big_midpoints(Kmi, 400, midpoints, which+4);
	sequence_midpoints(Kmi, which+4);
	double mmi = midpoints[which];

printf("Kmi=%g   mmi=%g\n", Kmi, mmi);

	if (mmi > mhi) return bisect(Kmi, Khi, which);
	if (mlo < mmi) { printf("Errrror !!!!!!!!!\n"); exit(1); }
	return bisect(Klo, Kmi, which);
}

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s maxn\n", argv[0]);
		exit(1);
	}
	int maxn = atoi(argv[1]);
	double guess = atof(argv[2]);

	double bad = bisect(guess-0.01, guess+0.01, maxn);

	printf("Its %d	%g\n", maxn, bad);

#ifdef PRINT_MIDPOINTS
#define NPTS 8701
	for (int i=0; i<NPTS; i++)
	{
		double x = ((double) i + 0.5) / ((double) NPTS);
		double K = 0.5 + 0.5*x;

		// find_midpoints(K);
		big_midpoints(K, 400, midpoints, maxn);
		sequence_midpoints(K, maxn);

		printf("%d	%g", i, K);
		for (int m=0; m< maxn; m++)
		{
			printf("	%g", midpoints[m]);
		}
		printf("\n");
	}
#endif

}
