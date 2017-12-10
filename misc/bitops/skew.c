
/*
 * skew.c
 *
 * Skew takago map
 * Dec 2017
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

int main (int argc, char* argv[])
{
	double K = atof(argv[1]);
	double w = atof(argv[2]);

#define NPTS 1803
	for (int i=0; i<NPTS; i++)
	{
		double x = ((double) i) / ((double) NPTS);
		double y = skew_takagi(x, K, w, 14);
		double z = skew_takagi(x, K, w, -1);
		double t = skew_takagi(x, K, w, 0);
		double s = skew_takagi(x, K, w, 1);
		double r = skew_takagi(x, K, w, 2);
		printf("%d	%g %g	%g	%g	%g	%g\n", i, x, y, z, t, s, r);
	}
}
