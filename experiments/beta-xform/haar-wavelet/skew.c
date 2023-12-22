/*
 * skew.c
 *
 * Skew Takagi map and skew Haar wavelet
 *
 * Dec 2017 Linas Vepstas
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double skew_triangle(double x, double K)
{
	x -= floor(x);
	if (x < 1.0/(2.0*K)) return  x * 2.0 * K;
	return (1.0 - x) / (1.0 - 1.0/(2.0*K));
}

double skew_square(double x, double K)
{
	x -= floor(x);
	if (x < 1.0/(2.0*K)) return  2.0 * K;
	return -1.0 / (1.0 - 1.0/(2.0*K));
}

double skew_takagi(double x, double K, double w, int cnt)
{
	if (cnt < 0)
		return skew_triangle(x,K);

	if (x < 1.0/(2.0*K)) 
	{
		double rescale = 2.0 * x * K;
		return rescale + w*skew_takagi(rescale, K, w, cnt-1);
	}
	else
	{
		double deno = 1.0 - 1.0/(2.0*K);
		return (1.0 - x) / deno +
			w*skew_takagi((x - 1.0/(2.0*K))/deno, K, w, cnt-1);
	}
}

double skew_haar(double x, double K, double w, int cnt)
{
	if (cnt < 0)
		return skew_square(x,K);

	if (x < 1.0/(2.0*K)) 
	{
		double rescale = 2.0 * x * K;
		return 2.0*K + w*skew_haar(rescale, K, w, cnt-1);
	}
	else
	{
		double deno = 1.0 - 1.0/(2.0*K);
		return - 1.0 / deno +
			w*skew_haar((x - 1.0/(2.0*K))/deno, K, w, cnt-1);
	}
}

int main (int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K w\n", argv[0]);
		exit (1);
	}
	double K = atof(argv[1]);
	double w = atof(argv[2]);

#define NPTS 1803
	for (int i=0; i<NPTS; i++)
	{
		double x = ((double) i) / ((double) NPTS);
#ifdef BASIC
		double y = skew_takagi(x, K, w, 14);
		double z = skew_takagi(x, K, w, -1);
		double t = skew_takagi(x, K, w, 0);
		double s = skew_takagi(x, K, w, 1);
		double r = skew_takagi(x, K, w, 2);
		printf("%d	%g %g	%g	%g	%g	%g\n", i, x, y, z, t, s, r);
#endif
#ifdef SELF_SIM
		double y = skew_takagi(x, K, w, 14);
		double lo = skew_takagi(0.5*x/K, K, w, 14);
		double hi = skew_takagi(0.5/K + x*(1.0 - 0.5/K), K, w, 14);
		printf("%d	%g %g	%g	%g\n", i, x, y, lo, hi);
#endif
#define BASIC_HAAR
#ifdef BASIC_HAAR
		double y = skew_haar(x, K, w, 14);
		double z = skew_haar(x, K, w, -1);
		double t = skew_haar(x, K, w, 0);
		double s = skew_haar(x, K, w, 1);
		double r = skew_haar(x, K, w, 2);
		printf("%d	%g %g	%g	%g	%g	%g\n", i, x, y, z, t, s, r);
#endif
#ifdef HAAR_SELF_SIM
		double y = skew_haar(x, K, w, 14);
		double lo = skew_haar(0.5*x/K, K, w, 14);
		double hi = skew_haar(0.5/K + x*(1.0 - 0.5/K), K, w, 14);
		printf("%d	%g %g	%g	%g\n", i, x, y, lo, hi);
#endif
	}
}
