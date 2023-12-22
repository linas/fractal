/*
 * psidist.c
 * Verify the distribution of the midpoints.
 * Generate the mid-points, and bin them. The result is as it should
 * be -- its nothing other than the staionary measure.
 *
 * February 2018
 */
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXN 400000
#define NOMAIN
#include "psi.c"
#include "psibig.c"

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K nbins\n", argv[0]);
		exit(1);
	}
	double K = atof(argv[1]);
	int nbins = atoi(argv[2]);
	printf("#\n# K=%g\n#\n", K);

	// find_midpoints(K);
	big_midpoints(K, 400, midpoints, MAXN);

	double s[nbins];
	for (int i=0; i< nbins; i++) s[i] = 0;

	for (int i=0; i<MAXN; i++)
	{
		double m = midpoints[i];
		// m /= K;
		m *= nbins;
		m = floor(m);
		int n = m;
		s[n] += 1.0;
	}
	for (int i=0; i< nbins; i++) s[i] *= nbins/((double) MAXN);

	// and again, high precision
	big_midpoints(K, 13200, midpoints, MAXN);

	double t[nbins];
	for (int i=0; i< nbins; i++) t[i] = 0;

	for (int i=0; i<MAXN; i++)
	{
		double m = midpoints[i];
		// m /= K;
		m *= nbins;
		m = floor(m);
		int n = m;
		t[n] += 1.0;
	}
	for (int i=0; i< nbins; i++) t[i] *= nbins/((double) MAXN);

	for (int i=0; i< nbins; i++)
	{
		double x = ((double) i + 0.5) / ((double) nbins);
		printf("%d	%g	%g	%g\n", i, x, s[i], t[i]);
	}
}
